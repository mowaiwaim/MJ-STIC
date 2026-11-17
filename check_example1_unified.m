function result = check_example1_unified()
% Examples 1 and 2 in NN_STIC_NEW. Evaluate all three passivity certificates.
% Also check the proposed threshold and candidate/switch arbitration on
% the reference states retained from the paper's two representative paths.
N = 5;
n = 3;
I = eye(N*n);
C1 = diag([1.7,1.8,1.9]);
C2 = diag([1.8,1.9,2.0]);
A1 = diag([10,5,12.5]);
A2 = diag([.3,.2,.3]);
B1 = [.1 .2 .3;.2 .4 .3;-.1 .3 .2];
B2 = [.1 .2 .2;.1 .3 .3;-.1 .3 .1];
E = kron(eye(N),[1 0 0;0 .9 0;1.1 0 .8;0 1.2 .7]);
D = kron(eye(N),[1.2 0 0;0 .9 1.1;0 0 1;1.2 0 0]);
Pi = kron(eye(N)-ones(N)/N,eye(n));
edges = [1 2;1 3;2 4;3 5;4 5];
gH = [.1 .2 .2;.2 .2 .3;.1 .2 .2;.2 .1 .2;.2 .1 .1];
gW = [.1 .1 .2;.1 .2 .1;.1 .2 .1;.2 .1 .2;.2 .1 .1];
psi = [.2 .3 .6;.8 .2 .4;.3 .4 .5];
a = [.1 .1 .2];
bd = [.3 .4 .5];
M = I;
Kc = zeros(size(I));
W = cell(3,1);
for q=1:3
    W{q} = kron(graph_matrix(N,edges,gW(:,q)),diag(psi(q,:)));
    Kc = Kc + a(q)*kron(graph_matrix(N,edges,gH(:,q)),diag(psi(q,:)));
    M = M - bd(q)*W{q};
end
rho = 1.32;
lambda = .01;
Hmax = 1;
nu = .002;
lowerWeight = exp(-nu*Hmax);
R = 1.30*I;
QfI = .5*max(eig(M*kron(eye(N),A2)));
pin = kron(diag([15 0 0 0 15]),eye(n));
T = M*(-kron(eye(N),C2)+Kc);
AI = T+T' + QfI*I + .5*M*kron(eye(N),B2*B2')*M + rho*.5*I;

bP = [15.2387 17.5570 15.8809;15.2915 17.9902 15.9517];
Z0 = {[189.9837 -95.4905 -187.3298;-259.5853 133.8342 259.6982; ...
       -171.6057 86.9009 170.4377;195.2924 -99.4824 -194.1457], ...
      [198.7424 -99.9351 -196.0024;-271.6738 140.0124 271.7423; ...
       -179.5495 90.9410 178.3383;204.3750 -104.1025 -203.1713]};
L0 = {zeros(4),[9.0771 -12.5602 -8.2427 9.4201; ...
      -12.5602 17.4091 11.4182 -13.0492; ...
      -8.2427 11.4182 7.5033 -8.5636;9.4201 -13.0492 -8.5636 9.7968]};
endpoint = zeros(2,3);
Lambda1 = .005*I;
inputStrictEndpoint = zeros(1,3);
for k=1:2
    Z = kron(eye(N),Z0{k});
    L = kron(eye(N),L0{k});
    AP = kron(eye(N),-2*C1+.5*A1+B1*B1') + rho*.5*I + 2*Kc - 2*pin;
    for q=1:3
        AP = AP+2*bP(k,q)*W{q};
    end
    % Settings are [mode, weight, clock decay, exit rate].
    clocks = [1 1 0 2;2 lowerWeight nu 8;2 1 nu 8];
    stateBlocks = {AP, AI};
    crossBlocks = {Pi, M*Pi};
    for j=1:3
        mode = clocks(j,1); p = clocks(j,2);
        X = p*(stateBlocks{mode}+(lambda-clocks(j,3))*M) ...
            + R + clocks(j,4)*(1-p)*M + 2*E'*L*E;
        U = p*crossBlocks{mode} - E'*Z + 2*E'*L*D;
        V = 2*D'*L*D - D'*Z - Z'*D;
        S = [X,U;U',V];
        endpoint(k,j) = max(eig((S+S')/2));
        if k==1
            % The input-strict test uses the ordinary gains, Z and clock.
            Sinput = [X,U;U',V+2*Lambda1];
            inputStrictEndpoint(j) = max(eig((Sinput+Sinput')/2));
        end
    end
end
delay = rho*.5+[.25,.50]-(1-.1)*exp(-lambda*.1)*1.30;
jump = (1-.2)^2*exp(nu*Hmax)-1;
assert(max(eig(M))<rho && min(eig(M))>0 && min(eig(L0{2}))>0)
assert(all(endpoint<0,'all') && all(delay<0) && jump<0)
assert(all(inputStrictEndpoint<0) && min(eig(Lambda1))>0)
disp('Maximum eigenvalues of the three passivity blocks:')
disp(array2table([endpoint(1,:);inputStrictEndpoint;endpoint(2,:)], ...
    'VariableNames',{'Available','ConstrainedLower','ConstrainedUpper'}, ...
    'RowNames',{'Ordinary','InputStrict','OutputStrict'}))

% Sampled energies and times from the two original simulation logs. An
% uncapped candidate beyond the old reporting horizon is evaluated here as
% a finite root and compared to the next switch on the same grid.
% Columns are time of reference, sampled Q, next switch and on-grid candidate.
logs{1} = [ ...
    0  1.421599  0.204  0.070;
    0.07  0.69761340735671  0.204  0.182;
    0.182  0.316532626920952  0.204  0.357;
    0.48  0.0341937310663538  0.542  0.892;
    0.575  0.054282969451236  0.581  0.934;
    1.111  0.00128905164811295  1.443  1.913;
    2.462  0.0159954680886248  2.686  2.963;
    2.922  0.00659091327914109  2.956  3.529;
];
logs{2} = [ ...
    0  1.421599  0.204  0.023;
    0.023  0.833222880317985  0.204  0.055;
    0.055  0.479388130238142  0.204  0.099;
    0.099  0.264077239300011  0.204  0.157;
    0.157  0.135056739200952  0.204  0.233;
    0.48  0.0142292427325699  0.542  0.624;
    0.575  0.00572736850224025  0.581  0.749;
    1.111  0.000853692043031045  1.443  1.347;
    1.347  0.000289487443657744  1.443  1.619;
    2.462  1.90434469491235e-08  2.686  3.055;
    2.922  1.31879242900636e-09  2.956  3.604;
];
schedules = [8 .05 .3;20 .01 10];
waitMax = zeros(2,1);
floorQ = zeros(2,1);
margins = zeros(2,1);
actual = cell(2,1);
dt = .001;
for k=1:2
    log = logs{k};
    ell = schedules(k,1); varpi = schedules(k,2); chi = schedules(k,3);
    floorQ(k) = exp(-chi*Hmax)/(expm1(ell*Hmax)-varpi);
    margins(k) = min(log(:,2))-floorQ(k);
    assert(margins(k)>0);
    delays = arrayfun(@(energy) fzero(@(d) max(energy,floorQ(k))* ...
        (expm1(ell*d)-varpi)-exp(-chi*d),[0,Hmax]),log(:,2));
    waitMax(k) = max(delays);
    onGrid = log(:,1)+dt*ceil(delays/dt);
    assert(max(abs(onGrid-log(:,4)))<1e-11);
    actual{k} = onGrid(onGrid<log(:,3));
end
assert(numel(actual{1})==2 && numel(actual{2})==5)
disp('Threshold, smallest sampled energy minus threshold, and largest waiting time:')
disp(table(floorQ,margins,waitMax))
result = struct('endpoint',endpoint,'inputStrictEndpoint',inputStrictEndpoint, ...
    'delay',delay,'jumpResidual',jump, ...
    'minimumEigenvalueOutputStrict',min(eig(L0{2})),'floorQ',floorQ, ...
    'minimumEnergyMargin',margins,'maximumWaiting',waitMax,'actual',{actual});
end

function g = graph_matrix(n,edges,weights)
g = zeros(n);
for k=1:size(edges,1)
    i=edges(k,1); j=edges(k,2);
    g(i,j)=weights(k); g(j,i)=weights(k);
end
g = g-diag(sum(g,2));
end
