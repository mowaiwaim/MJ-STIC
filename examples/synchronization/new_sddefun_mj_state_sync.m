function [driftX, driftS, sigmaX, sigmaS] = new_sddefun_mj_state_sync(t, X, Xd, s, sd, mode, Alpha)
% One-step drift/diffusion builder for the Markov jump hybrid system.
% mode = 1: continuous implicit leader-pinning PD active, no impulsive control.
% mode = 2: impulsive mode, no continuous PD.

persistent n N a Psi H W C A B bP bD MP MD MPbar Ppin IMD
if isempty(n)
    n = 3; N = 5;
    a = [0.1 0.1 0.2];

    Psi = cell(1,3);
    Psi{1} = diag([0.2 0.3 0.6]);
    Psi{2} = diag([0.8 0.2 0.4]);
    Psi{3} = diag([0.3 0.4 0.5]);

    H = cell(1,3);
    H{1} = [-0.3 0.1 0.2 0   0;
            0.1 -0.2 0   0.1 0;
            0.2 0   -0.4 0   0.2;
            0   0.1 0   -0.3 0.2;
            0   0   0.2 0.2 -0.4];
    H{2} = [-0.4 0.2 0.2 0   0;
            0.2 -0.4 0   0.2 0;
            0.2 0   -0.3 0   0.1;
            0   0.2 0   -0.3 0.1;
            0   0   0.1 0.1 -0.2];
    H{3} = [-0.5 0.2 0.3 0   0;
            0.2 -0.4 0   0.2 0;
            0.3 0   -0.5 0   0.2;
            0   0.2 0   -0.3 0.1;
            0   0   0.2 0.1 -0.3];

    W = cell(1,3);
    W{1} = [-0.2 0.1 0.1 0   0;
             0.1 -0.2 0   0.1 0;
             0.1 0   -0.3 0   0.2;
             0   0.1 0   -0.3 0.2;
             0   0   0.2 0.2 -0.4];
    W{2} = [-0.3 0.1 0.2 0   0;
             0.1 -0.3 0   0.2 0;
             0.2 0   -0.3 0   0.1;
             0   0.2 0   -0.3 0.1;
             0   0   0.1 0.1 -0.2];
    W{3} = [-0.3 0.2 0.1 0   0;
             0.2 -0.3 0   0.1 0;
             0.1 0   -0.3 0   0.2;
             0   0.1 0   -0.2 0.1;
             0   0   0.2 0.1 -0.3];

    C = cell(1,2);
    C{1} = diag([1.7 1.8 1.9]);
    C{2} = diag([1.8 1.9 2.0]);

    A = cell(1,2);
    A{1} = diag([10 5 12.5]);
    A{2} = diag([0.3 0.2 0.3]);

    B = cell(1,2);
    B{1} = [ 0.1 0.2 0.3;
             0.2 0.4 0.3;
            -0.1 0.3 0.2];
    B{2} = [ 0.1 0.2 0.2;
             0.1 0.3 0.3;
            -0.1 0.3 0.1];

    bD = [0.3 0.4 0.5];
    bP = [15.2915 17.9902 15.9517];

    MP = zeros(n*N);
    MD = zeros(n*N);
    for q = 1:3
        MP = MP + bP(q) * kron(W{q}, Psi{q});
        MD = MD + bD(q) * kron(W{q}, Psi{q});
    end
    DP    = diag([15 0 0 0 15]);
    Ppin  = kron(DP, eye(n));
    MPbar = MP - Ppin;
    IMD = eye(n*N) - MD;
end

Ck = C{mode};
Ak = A{mode};
Bk = B{mode};

fX  = f_act(X);
fS  = f_act(s);
gXd = g_act(Xd);
gSd = g_act(sd);

BaseX = -Ck * X + Ak * fX + Bk * gXd;
for q = 1:3
    BaseX = BaseX + a(q) * (Psi{q} * X) * H{q}.';
end
BaseX = BaseX + Alpha;

alpha_bar = mean(Alpha, 2);
BaseS = -Ck * s + Ak * fS + Bk * gSd + alpha_bar;

if mode == 1
    xvec    = reshape(X, [], 1);
    hatxvec = reshape(BaseX, [], 1);
    svec    = repmat(s, N, 1);
    uvec    = IMD \ (MPbar * xvec + MD * hatxvec + Ppin * svec);
    U       = reshape(uvec, n, N);
else
    U = zeros(n, N);
end

driftX = BaseX + U;
driftS = BaseS;

if mode == 1
    sigmaX = 0.5 * X + 0.5 * Xd;
    sigmaS = 0.5 * s + 0.5 * sd;
else
    sigmaX = 0.5 * X - 0.5 * Xd;
    sigmaS = 0.5 * s - 0.5 * sd;
end

if nargin < 1 %#ok<UNRCH>
    disp(t);
end
end

function y = f_act(x)
y = 0.25 * x - x.^3;
end

function y = g_act(x)
y = 0.5 * tanh(x);
end
