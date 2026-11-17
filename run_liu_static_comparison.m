function run_liu_static_comparison(sourceDir, outputDir)
% Reproduce Eq. (5) of Liu et al. (2026), doi:10.3390/math14111787,
% without attacks, with one timer per node during visits to mode 2.
if nargin < 1, sourceDir = fileparts(mfilename('fullpath')); end
if nargin < 2, outputDir = sourceDir; end
addpath(sourceDir);
p = bw_parameters();
pilot = p; pilot.samples = 256; pilot.seed = 512091;
paths = bw_paths(pilot);
ref = simulate_bounded_wait(pilot,paths,'stic');
target = numel(p.x0)*mean(ref.impulses);
assert(abs(target-48.26953125)<1e-10);
% Select gamma using only mean node-correction counts. Search by bracketing
% on log scale and keep the closest observed count (smaller scale breaks ties).
lo = -9; hi = 1;
scales = [10^lo;10^hi];
counts = zeros(size(scales));
for k = 1:2
    trial = simulate(pilot,paths,scales(k));
    counts(k) = mean(trial.corrections);
end
assert(counts(1)>target && counts(2)<target);
for k = 1:12
    mid = 0.5*(hi+lo);
    trial = simulate(pilot,paths,10^mid);
    scales(end+1,1) = 10^mid; %#ok<AGROW>
    counts(end+1,1) = mean(trial.corrections); %#ok<AGROW>
    if counts(end)>target, lo = mid; else, hi = mid; end
end
pilotTable = sortrows(table(scales,counts,repmat(target,size(counts)), ...
    'VariableNames',{'Scale','MeanNodeCorrections','Target'}),'Scale');
[~,choice] = min(abs(pilotTable.MeanNodeCorrections-target));
pilotTable.Selected = false(height(pilotTable),1); pilotTable.Selected(choice) = true;
scale = pilotTable.Scale(choice);
fprintf('Pilot target %.7f; selected scale %.12g; observed %.7f\n',target,scale,pilotTable.MeanNodeCorrections(choice));
writetable(pilotTable,fullfile(outputDir,'liu_pilot.csv'));
test = simulate(p,bw_paths(p),scale);
data = table((1:p.samples)',test.energy(end,:)',trapz(test.time,test.energy,1)', ...
    test.impulses,test.corrections,NaN(p.samples,1),test.samples,test.capped, ...
    test.byRule,'VariableNames',{'Trial','TerminalSquaredError','IntegralSquaredError', ...
    'Impulses','NodeCorrections','FullErrorSnapshots','ScalarErrorSamples', ...
    'UpperBoundSelections','RuleSelections'});
% Impulses counts distinct times of applied corrections (ties grouped).
% ScalarErrorSamples counts one pre-jump error per actuated node and one
% local sample per node at each mode-2 entry, including initialization.
writetable(data,fullfile(outputDir,'liu_paths.csv'));
summary=table("Liu et al. (2026), static",scale,mean(data.Impulses), ...
    mean(data.NodeCorrections),NaN,mean(data.ScalarErrorSamples),mean(data.IntegralSquaredError), ...
    mean(data.TerminalSquaredError),mean(data.UpperBoundSelections),mean(data.RuleSelections), ...
    'VariableNames',{'Method','Scale','Impulses','NodeCorrections','FullErrorSnapshots','ScalarErrorSamples', ...
    'IntegratedMSE','TerminalMSE','UpperBoundSelections','RuleSelections'});
writetable(summary,fullfile(outputDir,'liu_summary.csv')); disp(summary);
curve=table(test.time,mean(test.energy,2),'VariableNames',{'Time','LiuStatic'});
writetable(curve,fullfile(outputDir,'liu_trajectories.csv'));
end

function result = simulate(p,paths,scale)
n = round(p.T/p.dt);
tgrid = (0:n)'*p.dt;
N = numel(p.x0);
I = eye(N);
MP = p.M\(-p.C(1)*I+p.a*p.Hgraph+p.bp*p.Wgraph-p.pin);
AP = p.M\(p.A(1)*I); BP = p.M\(p.B(1)*I);
MI = -p.C(2)*I+p.a*p.Hgraph;
energy = zeros(n+1,p.samples);
impulses = zeros(p.samples,1); corrections = impulses;
samples = impulses; capped = impulses; byRule = impulses;
for trial=1:p.samples
    x=p.x0; X=zeros(N,n+1); X(:,1)=x;
    times=[paths.switches{trial},inf]; mode=2; ns=1;
    waits=localWaits(x,0,scale,p); deadline=waits;
    samples(trial)=N;
    capped(trial)=sum(waits==p.H); byRule(trial)=sum(waits<p.H);
    for k=1:n
        now=tgrid(k); finish=tgrid(k+1); db=paths.dB(k,trial);
        while now<finish
            due=min(deadline);
            tNext=min([finish,times(ns),due]); h=tNext-now;
            if h>0
                lagIndex=max(1,1+floor((now-p.delay)/p.dt+1e-11));
                lag=X(:,lagIndex);
                if mode==1
                    drift=MP*x+AP*tanh(x)+BP*(p.gSlope*tanh(lag));
                else
                    drift=MI*x+p.A(2)*tanh(x)+p.B(2)*p.gSlope*tanh(lag);
                end
                x=x+h*drift+p.noise*x*(h/p.dt*db);
            end
            if tNext==times(ns)
                ns=ns+1; mode=3-mode;
                deadline(:)=inf;
                if mode==2
                    waits=localWaits(x,tNext,scale,p);
                    deadline=tNext+waits; samples(trial)=samples(trial)+N;
                    capped(trial)=capped(trial)+sum(waits==p.H);
                    byRule(trial)=byRule(trial)+sum(waits<p.H);
                end
            elseif tNext==due
                dueNodes=(deadline==due);
                % Eq. (5) samples the corresponding node before its jump.
                waits=localWaits(x(dueNodes),tNext,scale,p);
                deadline(dueNodes)=tNext+waits;
                x(dueNodes)=(1+p.mu)*x(dueNodes);
                impulses(trial)=impulses(trial)+1;
                corrections(trial)=corrections(trial)+sum(dueNodes);
                samples(trial)=samples(trial)+sum(dueNodes);
                capped(trial)=capped(trial)+sum(waits==p.H);
                byRule(trial)=byRule(trial)+sum(waits<p.H);
            end
            now=tNext;
        end
        X(:,k+1)=x;
    end
    energy(:,trial)=sum(X.^2,1)';
end
assert(all(samples==capped+byRule));
result=struct('energy',energy,'time',tgrid,'impulses',impulses, ...
    'corrections',corrections,'samples',samples,'capped',capped,'byRule',byRule);
end

function waits=localWaits(x,t,scale,p)
% Multiplying Phi by exp(-alpha2*(t+d)) preserves its sign and makes
% every term increasing for d>0 because alpha1>alpha2; bisection finds
% the *first* crossing, and the value at zero is never evaluated.
kappa=.05; alpha1=5; alpha2=.15;
g1=1.7*scale; g2=1.5*scale;
energy=x.^2;
waits=repmat(p.H,size(x));
positive = margin(energy,t,p.H,kappa,alpha1,alpha2,g1,g2)>0;
if any(positive)
    E=energy(positive);
    lo = repmat(log1p(kappa)/alpha1,size(E)); hi = repmat(p.H,size(E));
    for k=1:38
        mid=0.5*(lo+hi);
        hit=margin(E,t,mid,kappa,alpha1,alpha2,g1,g2)>0;
        hi(hit)=mid(hit); lo(~hit)=mid(~hit);
    end
    waits(positive)=0.5*(hi+lo);
end
end

function value=margin(E,t,d,kappa,alpha1,alpha2,g1,g2)
value=E.*(expm1(alpha1*d)-kappa)-g1*exp(alpha2*(t+d))-g2./(t+d).^2;
end
