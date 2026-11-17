function run_shang_comparison(applyCap)
% Independent pulse-count calibration and paired test of Shang & Zhang (2025).
% Required on the path: bw_parameters, bw_paths, simulate_bounded_wait.
% Source: AIMS Mathematics 10, 20368-20384, doi:10.3934/math.2025910.
% We use the two-sample rule (18), the waiting function in Example 2 and
% the squared error energy of Example 3 in the present manuscript.
if nargin==0, applyCap=true; end
if applyCap, folder='outputs'; else, folder='uncapped'; end
outdir=fullfile(fileparts(mfilename('fullpath')),folder);
if ~isfolder(outdir), mkdir(outdir); end
p = bw_parameters();
pilot = p;
pilot.seed = 512091;
pilot.samples = 256;
trainPaths = bw_paths(pilot);
trainMJ = simulate_bounded_wait(pilot,trainPaths,'stic',NaN);
target = mean(trainMJ.impulses);
% Search only pulse counts. Theta=2, max(pi)=0.25, xi=1/4; this
% range respects theta*max(pi)<ell and ell+2*0.01*max(pi)<-log(xi).
lo = 0.501;
hi = 1.35;
trials = zeros(14,2);
for j = 1:size(trials,1)
    if j==1
        ell=lo;
    elseif j==2
        ell=hi;
    else
        ell=(lo+hi)/2;
    end
    sim=simulate_shang(pilot,trainPaths,ell,applyCap);
    n=mean(sim.impulses);
    trials(j,:)=[ell,n];
    if j==2
        assert(trials(1,2)>=target && n<=target,'Pilot budget outside preset design range');
    elseif j>2
        if n>target, lo=ell; else, hi=ell; end
    end
end
[~,best]=min(abs(trials(:,2)-target));
ell=trials(best,1);
scores=array2table(trials,'VariableNames',{'Ell','PilotPulses'});
scores.MJTarget=repmat(target,height(scores),1);
scores.Selected=(1:height(scores))'==best;
writetable(scores,fullfile(outdir,'shang_pilot.csv'));

testPaths=bw_paths(p);
run=simulate_shang(p,testPaths,ell,applyCap);
out=table((1:p.samples)',run.terminal,run.integral,run.impulses,numel(p.x0)*run.impulses, ...
    run.snapshots,run.references,run.probes,run.immediate,run.cancelProbe,run.cancelPulse,run.capped,run.cappedExecutions, ...
    'VariableNames',{'Trial','TerminalSquaredError','IntegralSquaredError','Impulses','NodeCorrections', ...
    'FullErrorSnapshots','ReferenceStarts','IntermediateSamples','CoincidentSamples','CancelledIntermediate','CancelledPulse','CappedPlans','CappedExecutions'});
writetable(out,fullfile(outdir,'shang_paths.csv'));
summary=array2table(mean(out{:,2:end},1),'VariableNames',out.Properties.VariableNames(2:end));
summary.CapSeconds=inf;
if applyCap, summary.CapSeconds=p.H; end
summary.CapShareOfPlans=sum(run.capped)/sum(run.probes);
summary.CapShareOfExecutions=sum(run.cappedExecutions)/sum(run.impulses);
summary.MinPlannedWait=min(run.minPlanned); summary.MaxPlannedWait=max(run.maxPlanned);
summary.MinExecutedWait=min(run.minExecuted); summary.MaxExecutedWait=max(run.maxExecuted);
summary.MaxUntruncatedPlannedWait=max(run.maxRaw);
summary.Ell=ell; summary.PilotTarget=target; summary.PilotPulses=trials(best,2);
writetable(summary,fullfile(outdir,'shang_summary.csv'));
curves=table(run.time,mean(run.energy,2),'VariableNames',{'Time','MeanSquaredError'});
writetable(curves,fullfile(outdir,'shang_trajectories.csv'));
disp(summary);
end

function out=simulate_shang(p,paths,ell,applyCap)
n=round(p.T/p.dt);
grid=(0:n)'*p.dt;
N=numel(p.x0);
I=eye(N);
MP=p.M\(-p.C(1)*I+p.a*p.Hgraph+p.bp*p.Wgraph-p.pin);
AP=p.M\(p.A(1)*I);
BP=p.M\(p.B(1)*I);
MI=-p.C(2)*I+p.a*p.Hgraph;
energy=zeros(n+1,p.samples);
impulses=zeros(p.samples,1); references=impulses; probes=impulses;
cancelProbe=impulses; cancelPulse=impulses; immediate=impulses;
capped=impulses; cappedExecutions=impulses; maxPlanned=impulses; maxRaw=impulses;
minPlanned=inf(p.samples,1); minExecuted=minPlanned; maxExecuted=impulses;
for trial=1:p.samples
    X=zeros(N,n+1); x=p.x0; X(:,1)=x;
    mode=2; ns=1; switches=[paths.switches{trial},inf];
    refTime=0; capActive=false; refEnergy=x'*p.M*x;
    nextProbe=waiting(x); deadline=inf;
    references(trial)=1;
    for k=1:n
        now=grid(k);
        finish=grid(k+1);
        while now<finish
            tNext=min([finish,switches(ns),nextProbe,deadline]);
            h=tNext-now;
            if h>0
                lagIndex=max(1,1+floor((now-p.delay)/p.dt+1e-11));
                xLag=X(:,lagIndex);
                if mode==1
                    drift=MP*x+AP*tanh(x)+BP*(p.gSlope*tanh(xLag));
                else
                    drift=MI*x+p.A(2)*tanh(x)+p.B(2)*p.gSlope*tanh(xLag);
                end
                x=x+h*drift+p.noise*x*(h/p.dt*paths.dB(k,trial));
            end
            if tNext==switches(ns)
                if mode==2
                    if isfinite(nextProbe)
                        cancelProbe(trial)=cancelProbe(trial)+1;
                    else
                        cancelPulse(trial)=cancelPulse(trial)+1;
                    end
                end
                ns=ns+1; refTime=tNext; capActive=false;
                mode=3-mode;
                deadline=inf; nextProbe=inf;
                if mode==2
                    refEnergy=x'*p.M*x;
                    nextProbe=tNext+waiting(x);
                    references(trial)=references(trial)+1;
                end
            else
                if tNext==nextProbe
                    probes(trial)=probes(trial)+1;
                    current=x'*p.M*x;
                    assert(refEnergy>0 && current>0,'Undefined energy ratio at an observed sample');
                    remaining=max(0,(ell+log(refEnergy/current))/2.01);
                    rawDeadline=tNext+remaining;
                    maxRaw(trial)=max(maxRaw(trial),rawDeadline-refTime);
                    capActive=applyCap && rawDeadline>refTime+p.H;
                    if capActive
                        deadline=refTime+p.H;
                        capped(trial)=capped(trial)+1;
                    else
                        deadline=rawDeadline;
                    end
                    assert(deadline>=tNext);
                    minPlanned(trial)=min(minPlanned(trial),deadline-refTime);
                    maxPlanned(trial)=max(maxPlanned(trial),deadline-refTime);
                    nextProbe=inf;
                    if deadline==tNext, immediate(trial)=immediate(trial)+1; end
                end
                if tNext==deadline
                    minExecuted(trial)=min(minExecuted(trial),tNext-refTime);
                    maxExecuted(trial)=max(maxExecuted(trial),tNext-refTime);
                    cappedExecutions(trial)=cappedExecutions(trial)+capActive;
                    refTime=tNext; capActive=false;
                    x=(1+p.mu)*x;
                    impulses(trial)=impulses(trial)+1;
                    references(trial)=references(trial)+1;
                    refEnergy=x'*p.M*x;
                    nextProbe=tNext+waiting(x); deadline=inf;
                end
            end
            now=tNext;
        end
        X(:,k+1)=x;
    end
    energy(:,trial)=sum(X.^2,1)';
end
assert(min(minExecuted)>=.15-1e-9);
if applyCap, assert(max(maxPlanned)<=p.H+1e-9); end
assert(all(probes>=impulses) && all(references>=impulses));
out=struct('time',grid,'energy',energy,'terminal',energy(end,:)', ...
    'integral',trapz(grid,energy,1)','impulses',impulses,'references',references, ...
    'probes',probes,'snapshots',references+probes-immediate,'cancelProbe',cancelProbe, ...
    'cancelPulse',cancelPulse,'immediate',immediate,'capped',capped, ...
    'cappedExecutions',cappedExecutions,'maxPlanned',maxPlanned,'minPlanned',minPlanned, ...
    'minExecuted',minExecuted,'maxExecuted',maxExecuted,'maxRaw',maxRaw);
end

function sleepTime=waiting(x)
sleepTime=0.15+0.1*exp(-100*norm(x));
end

