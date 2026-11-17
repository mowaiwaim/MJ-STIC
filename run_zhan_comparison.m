function run_zhan_comparison()
% Reuse Example 3's network, paths and solver. Only the mode-2 timer changes.
% Zhan et al., IEEE/CAA J. Autom. Sinica 12 (2025), 264-266, Eq. (5).
% Run with bw_parameters, bw_paths, and simulate_bounded_wait on the path.
p = bw_parameters();
outDir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~isfolder(outDir), mkdir(outDir); end
rule.dwell = 0.1;
rule.growth = 2;
rule.rate = 0.01;
pilot = p;
pilot.samples = 256;
pilot.seed = 512091;
train = bw_paths(pilot);
trainBase = simulate_bounded_wait(pilot,train,'stic');
target = mean(trainBase.impulses);
assert(abs(target-16.08984375)<1e-12);
% Fit gamma solely to the number of executed pilot pulses. The endpoints
% meet gamma > growth*dwell and gamma + rate*dwell < -log((1+mu)^2).
lo = 0.21;
hi = 1.30;
rows = zeros(0,2);
for gamma = [lo,hi]
    tmp = simulate_rule(pilot,train,rule,gamma);
    rows(end+1,:) = [gamma,mean(tmp.impulses)]; %#ok<AGROW>
end
assert(rows(1,2)>=target && rows(2,2)<=target);
for k=1:13
    gamma=(lo+hi)/2;
    tmp=simulate_rule(pilot,train,rule,gamma);
    pulseMean=mean(tmp.impulses);
    rows(end+1,:)=[gamma,pulseMean]; %#ok<AGROW>
    if pulseMean==target, break; end
    if pulseMean>target, lo=gamma; else, hi=gamma; end
end
[~,selected]=min(abs(rows(:,2)-target));
gamma=rows(selected,1);
fit = array2table(rows,'VariableNames',{'Gamma','MeanPilotImpulses'});
fit.Target = repmat(target,size(rows,1),1);
fit.Selected = (1:height(fit))'==selected;
writetable(fit,fullfile(outDir,'zhan_pilot.csv'));
paths = bw_paths(p);
res = simulate_rule(p,paths,rule,gamma);
assert(all(res.maxAge<=p.H+1e-12));
assert(min(res.waits)<p.H-1e-3 && max(res.waits)>min(res.waits)+1e-3);
Trial=(1:p.samples)'; TerminalSquaredError=res.energy(end,:)';
IntegralSquaredError=trapz(res.time,res.energy,1)'; Impulses=res.impulses;
NodeCorrections=numel(p.x0)*Impulses;
FullErrorSnapshots=res.references+res.probes;
References=res.references; ProbeSamples=res.probes;
CanceledTimers=res.cancellations; UsedCap=res.capped;
perPath=table(Trial,TerminalSquaredError,IntegralSquaredError,Impulses, ...
    NodeCorrections,FullErrorSnapshots,References,ProbeSamples,CanceledTimers,UsedCap);
writetable(perPath,fullfile(outDir,'zhan_paths.csv'));
meanRow=varfun(@mean,perPath(:,2:end));
meanRow.Gamma=gamma; meanRow.Delta=rule.dwell;
meanRow.Lambda=rule.growth; meanRow.Lambda0=rule.rate;
meanRow.LowerWait=min(res.waits); meanRow.UpperWait=max(res.waits);
meanRow.UncappedFraction=mean(res.waits<p.H-1e-12);
writetable(meanRow,fullfile(outDir,'zhan_summary.csv'));
curves=table(res.time,mean(res.energy,2),'VariableNames',{'Time','MeanSquaredError'});
writetable(curves,fullfile(outDir,'zhan_trajectories.csv'));
disp(meanRow)
end

function out = simulate_rule(p,paths,rule,gamma)
% Same drift, EM stepping, delay discretization, off-grid Brownian interpolation
% and switch priority as simulate_bounded_wait.m. The additional probe at
% reference+dwell observes the full error without executing a control jump.
n=round(p.T/p.dt); tgrid=(0:n)'*p.dt; N=numel(p.x0);
assert(isequal(size(paths.dB),[n,p.samples]) && p.delay>=p.dt);
assert(rule.dwell>0 && rule.dwell<p.H && gamma>rule.growth*rule.dwell);
assert(gamma+rule.rate*rule.dwell < -log((1+p.mu)^2));
denom=rule.growth+rule.rate; I=eye(N);
MP=p.M\(-p.C(1)*I+p.a*p.Hgraph+p.bp*p.Wgraph-p.pin);
AP=p.M\(p.A(1)*I); BP=p.M\(p.B(1)*I); MI=-p.C(2)*I+p.a*p.Hgraph;
out.time=tgrid;out.energy=zeros(n+1,p.samples);
out.impulses=zeros(p.samples,1);out.references=ones(p.samples,1);
out.probes=zeros(p.samples,1);out.cancellations=zeros(p.samples,1);
out.capped=zeros(p.samples,1);out.maxAge=zeros(p.samples,1);
% Keep an actual multiset of announced waiting times for the held-out audit.
waitCells=cell(p.samples,1);
for trial=1:p.samples
    X=zeros(N,n+1); x=p.x0; X(:,1)=x; mode=2;
    switches=[paths.switches{trial},inf]; ns=1;
    reference=0; qRef=x'*p.M*x;
    % Every mode-2 entry uses the initial convention bar(zeta_0)=zeta_0:
    % both values in the logarithm are the entry sample.
    entryWait=min(p.H,rule.dwell+gamma/denom);
    deadline=entryWait; probe=inf; waits=entryWait;
    out.capped(trial)=double(entryWait==p.H);
    for k=1:n
        now=tgrid(k); finish=tgrid(k+1);
        while now<finish
            tNext=min([finish,switches(ns),deadline,probe]); step=tNext-now;
            if step>0
                lag=max(1,1+floor((now-p.delay)/p.dt+1e-11)); xLag=X(:,lag);
                if mode==1
                    drift=MP*x+AP*tanh(x)+BP*(p.gSlope*tanh(xLag));
                else
                    drift=MI*x+p.A(2)*tanh(x)+p.B(2)*p.gSlope*tanh(xLag);
                end
                x=x+step*drift+p.noise*x*(step/p.dt*paths.dB(k,trial));
            end
            if tNext==switches(ns)
                if mode==2
                    out.cancellations(trial)=out.cancellations(trial)+1;
                    out.maxAge(trial)=max(out.maxAge(trial),tNext-reference);
                end
                mode=3-mode; ns=ns+1; probe=inf; deadline=inf;
                if mode==2
                    reference=tNext; qRef=x'*p.M*x;
                    deadline=tNext+entryWait;
                    out.references(trial)=out.references(trial)+1;
                    out.capped(trial)=out.capped(trial)+(entryWait==p.H);
                    waits(end+1,1)=entryWait; %#ok<AGROW>
                end
            else
                if tNext==probe
                    qProbe=x'*p.M*x;
                    assert(qRef>0 && qProbe>0,'Logarithm needs nonzero sampled energies.');
                    predicted=(gamma+log(qRef/qProbe))/denom;
                    % The infimum in (5) equals the end of the dwell whenever
                    % the predicted extension is nonpositive.
                    d=min(p.H,rule.dwell+max(0,predicted));
                    deadline=max(tNext,reference+d); probe=inf;
                    out.probes(trial)=out.probes(trial)+1;
                    out.capped(trial)=out.capped(trial)+(d==p.H);
                    waits(end+1,1)=d; %#ok<AGROW>
                end
                if tNext==deadline
                    out.maxAge(trial)=max(out.maxAge(trial),tNext-reference);
                    x=(1+p.mu)*x;
                    out.impulses(trial)=out.impulses(trial)+1;
                    reference=tNext; qRef=x'*p.M*x;
                    out.references(trial)=out.references(trial)+1;
                    probe=tNext+rule.dwell; deadline=inf;
                end
            end
            now=tNext;
        end
        X(:,k+1)=x;
    end
    if mode==2, out.maxAge(trial)=max(out.maxAge(trial),p.T-reference); end
    out.energy(:,trial)=sum(X.^2,1)';
    waitCells{trial}=waits;
end
out.waits=vertcat(waitCells{:});
end

