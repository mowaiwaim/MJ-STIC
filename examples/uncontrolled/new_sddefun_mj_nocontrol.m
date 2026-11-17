function [driftX, driftS, sigmaX, sigmaS] = new_sddefun_mj_nocontrol(t, X, Xd, s, sd, mode, Alpha)
% One-step drift/diffusion builder for the Markov jump leader-follower system
% WITHOUT control.
% mode = 1 or 2 only changes the mode-dependent matrices and diffusion law.

persistent a Psi H C A B
if isempty(a)
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

driftX = BaseX;
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
