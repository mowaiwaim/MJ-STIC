function out = simulate_bounded_wait(p,paths,method,interval)
% Event-aligned Euler-Maruyama on shared CTMC and Brownian paths.
% At an off-grid event the fixed Brownian path is linearly interpolated.
% The delayed state is the stored, left-grid observation of t-delay.
if nargin < 4, interval = NaN; end
method = validatestring(method,{'stic','zero_floor','none','periodic','lu'});
periodic = any(strcmp(method,{'periodic','lu'}));
if periodic, assert(isfinite(interval) && interval>0); end
n = round(p.T/p.dt);
N = numel(p.x0);
nt = n+1;
assert(isequal(size(paths.dB),[n p.samples]) && p.delay >= p.dt);
tgrid = (0:n)'*p.dt;
energy = zeros(nt,p.samples);
impulses = zeros(p.samples,1);
cancellations = impulses;
references = impulses;
maxAge = impulses;
switchCount = cellfun(@numel,paths.switches);
firstEvents = zeros(0,3);
firstCandidates = zeros(0,3);
I = eye(N);
MP = p.M\(-p.C(1)*I+p.a*p.Hgraph+p.bp*p.Wgraph-p.pin);
AP = p.M\(p.A(1)*I);
BP = p.M\(p.B(1)*I);
MI = -p.C(2)*I+p.a*p.Hgraph;
for trial=1:p.samples
    X = zeros(N,nt);
    x=p.x0; X(:,1)=x;
    mode=2;
    switches=[paths.switches{trial},inf];
    ns=1;
    deadline = inf; reference = 0;
    if ~strcmp(method,'none')
        deadline = waiting(x,method,interval,p);
        references(trial)=1;
        if trial==1, firstCandidates(end+1,:)=[0,deadline,x'*p.M*x]; end
    end
    for k=1:n
        now=tgrid(k);
        finish=tgrid(k+1);
        while now < finish
            tNext=min([finish,switches(ns),deadline]);
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
                if mode==2 && ~strcmp(method,'none')
                    cancellations(trial)=cancellations(trial)+1;
                    maxAge(trial)=max(maxAge(trial),tNext-reference);
                end
                ns=ns+1;
                mode=3-mode;
                deadline=inf;
                if mode==2 && ~strcmp(method,'none')
                    reference=tNext;
                    deadline=tNext+waiting(x,method,interval,p);
                    references(trial)=references(trial)+1;
                    if trial==1, firstCandidates(end+1,:)=[reference,deadline,x'*p.M*x]; end %#ok<AGROW>
                end
            elseif tNext==deadline
                before=sum(x.^2);
                if strcmp(method,'lu')
                    [~,node]=max(abs(x)); % Ties go to the lowest index.
                    x(node)=(1+p.mu)*x(node);
                else
                    x=(1+p.mu)*x;
                end
                impulses(trial)=impulses(trial)+1;
                maxAge(trial)=max(maxAge(trial),tNext-reference);
                if trial==1, firstEvents(end+1,:)=[tNext,before,sum(x.^2)]; end %#ok<AGROW>
                reference=tNext;
                deadline=tNext+waiting(x,method,interval,p);
                references(trial)=references(trial)+1;
                if trial==1, firstCandidates(end+1,:)=[reference,deadline,x'*p.M*x]; end %#ok<AGROW>
            end
            now=tNext;
        end
        X(:,k+1)=x;
    end
    if mode==2 && ~strcmp(method,'none'), maxAge(trial)=max(maxAge(trial),p.T-reference); end
    energy(:,trial)=sum(X.^2,1)';
end
out.time=tgrid;
out.energy=energy;
out.terminal=energy(end,:)';
out.integral=trapz(tgrid,energy,1)';
out.impulses=impulses;
out.cancellations=cancellations;
out.references=references;
out.maxAge=maxAge;
out.switchCount=switchCount;
out.firstEvents=firstEvents;
out.firstCandidates=firstCandidates;
end

function d=waiting(x,method,interval,p)
if strcmp(method,'lu') || strcmp(method,'periodic'), d=interval; return; end
Q=x'*p.M*x;
if strcmp(method,'stic')
    sampledEnergy=max(Q,p.q);
    if sampledEnergy==p.q
        d=p.H; % Exact root at the positive floor, fixed at design time.
        return
    end
    hi=p.H;
else
    sampledEnergy=Q; % Zero-floor ablation only.
    if sampledEnergy==0, d=inf; return; end
    hi=log1p(p.varpi+1/sampledEnergy)/p.ell;
end
lo=log1p(p.varpi)/p.ell;
for i=1:42
    mid=(lo+hi)/2;
    if sampledEnergy*(expm1(p.ell*mid)-p.varpi)>exp(-p.chi*mid), hi=mid; else, lo=mid; end
end
d=(lo+hi)/2;
end
