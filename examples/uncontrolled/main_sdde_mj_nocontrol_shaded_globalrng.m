function main_sdde_mj_nocontrol_shaded_globalrng
% Markov jump SDDE leader-follower network WITHOUT control
% alpha_i(t)=0, continuous PD off, impulsive control off.
% This version keeps the OLD random-number behavior:
%   - use rng(2,'twister') once, matching the controlled example
%   - use rand for Markov jumps
%   - use randn for Wiener increments
% so Markov jumps and Wiener noise share the same global stream.
%
% Diagnostic figures only (not the paper Monte Carlo figure):
%   1) ||x_i(t)||_2 and ||s(t)||_2 with mode-2 shading
%   2) single-path leader/follower squared disagreement D_L and D_F
% The publication open/closed comparison is generated only by
% export_my8paper_figures.m, which uses paired common random numbers.

clc; close all;
clear new_sddefun_mj_nocontrol;
rng(2,'twister');

%% basic settings
 t0 = 0;
 tfinal = 10;
 dt = 1e-3;
 T = t0:dt:tfinal;
 Nt = numel(T);

 n = 3;
 N = 5;

 Shist = [0.3; -0.2; 0.4];
 Xhist = [ 0.02 -0.03  0.04 -0.05  0.06;
           0.01 -0.01  0.01 -0.02  0.02;
           0.05 -0.07  0.10 -0.13  0.16];
 X = zeros(n, N, Nt);
 S = zeros(n, Nt);
 X(:,:,1) = Xhist;
 S(:,1)   = Shist;

 % Same generator as the controlled and passivity simulations.
 Gamma = [-2 2; 8 -8];
 r0 = 2;
 r  = simulate_markov_chain_grid(Gamma, T, r0);
 fprintf('Markov path: %d switches, mode-2 fraction = %.4f\n', ...
     nnz(diff(r)),mean(r==2));

%% Euler-Maruyama time stepping
for k = 1:Nt-1
    tk = T(k);
    hk = new_delay(tk, []);
    tdk = tk - hk;

    Xd = get_delayed_followers(tdk, T, X, Xhist);
    Sd = get_delayed_leader(tdk, T, S, Shist);

    Alpha = alpha_input_zero(tk, N);
    [driftX, driftS, sigmaX, sigmaS] = new_sddefun_mj_nocontrol(tk, X(:,:,k), Xd, S(:,k), Sd, r(k), Alpha);

    % OLD style: Wiener increment from the same global RNG
    dW = sqrt(dt) * randn;

    X(:,:,k+1) = X(:,:,k) + driftX * dt + sigmaX * dW;
    S(:,k+1)   = S(:,k)   + driftS * dt + sigmaS * dW;
end

%% post-processing
hx = squeeze(sqrt(sum(X.^2,1)));
hs = sqrt(sum(S.^2,1));
E  = X - reshape(S, n, 1, Nt);
eh = squeeze(sqrt(sum(E.^2,1)));
Xbar = mean(X,2);
% D_L: node-averaged squared leader-tracking error for this sample path.
DL = squeeze(mean(sum(E.^2,1),2));
% D_F: node-averaged squared follower dispersion about their own mean.
DF = squeeze(mean(sum((X-Xbar).^2,1),2));
tail = T >= 0.9*tfinal;

fprintf('Open-loop D_L(0) = %.8e, D_L(T) = %.8e, tail mean = %.8e\n', ...
    DL(1), DL(end), mean(DL(tail)));
fprintf('Open-loop D_F(0) = %.8e, D_F(T) = %.8e, tail mean = %.8e\n', ...
    DF(1), DF(end), mean(DF(tail)));
fprintf('Leader norm: initial = %.8e, min = %.8e, max = %.8e, final = %.8e\n', ...
    hs(1),min(hs),max(hs),hs(end));
mode2_intervals = get_mode2_intervals(T, r);

%% plots
colors = [0.12 0.56 1.00;
          1.00 0.84 0.00;
          0.90 0.25 0.25;
          0.60 0.20 0.98;
          0.18 0.55 0.34];

figure(1); clf; hold on; grid on;
for i = 1:N
    plot(T, hx(i,:), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
end
plot(T, hs, 'k--', 'LineWidth', 1.6);
xlabel('time t(s)'); %ylabel('norm value');
%title('Open-loop state norms: follower states and leader state');
leg = cell(1, N+1);
for i = 1:N
    leg{i} = sprintf('$\\|x_%d(t)\\|_2$', i);
end
leg{N+1} = '$\|s(t)\|_2$';
legend(leg, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');
add_mode2_shading(ax(gca), mode2_intervals);
axx = gca;
axx.Box = 'on';
axx.XColor = 'k';
axx.YColor = 'k';
axx.LineWidth = 1.0;
hold off;

figure(2); clf; hold on; grid on;
semilogy(T, max(DL,realmin), 'b-', 'LineWidth', 1.4);
semilogy(T, max(DF,realmin), 'r--', 'LineWidth', 1.4);
xlabel('time t(s)');
ylabel('node-averaged squared disagreement');
legend({'$D_L(t)$','$D_F(t)$'}, 'Interpreter', 'latex', 'Location', 'northeast');
title('Uncontrolled network: persistent disagreement');
hold off;

%figure(2); clf; hold on; grid on;
%for i = 1:N
 %   plot(T, eh(i,:), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
%end
%xlabel('time t(s)'); ylabel('error norm');
%title('Open-loop synchronization errors');
%leg2 = cell(1, N);
%for i = 1:N
 %   leg2{i} = sprintf('$\\|x_%d(t)-s(t)\\|_2$', i);
%end
%legend(leg2, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
%add_mode2_shading(ax(gca), mode2_intervals);
%hold off;

%figure(3); clf;
%stairs(T, r, 'LineWidth', 1.2);
%grid on;
%xlabel('t'); ylabel('Markovian switching signal');
%ylim([0 2.5]); yticks([1 2]);
%title('Markovian switching signal r(t)');
end

function a = ax(a), end

function add_mode2_shading(ax, mode2_intervals)
yl = ylim(ax);
for k = 1:size(mode2_intervals,1)
    x1 = mode2_intervals(k,1);
    x2 = mode2_intervals(k,2);
    p = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(p,'bottom');
    xline(ax, x1, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, x2, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
ylim(ax, yl);
end

function intervals = get_mode2_intervals(T, r)
mask = (r == 2);
starts = find(diff([0 mask]) == 1);
stops  = find(diff([mask 0]) == -1);
intervals = zeros(numel(starts), 2);
for k = 1:numel(starts)
    i1 = starts(k); i2 = stops(k);
    t1 = T(i1);
    if i2 < numel(T)
        t2 = T(i2+1);
    else
        t2 = T(i2);
    end
    intervals(k,:) = [t1, t2];
end
end

function r = simulate_markov_chain_grid(Q, T, r0)
dt = T(2) - T(1);
P  = expm(Q * dt);
P  = max(P, 0);
P  = P ./ sum(P, 2);
Nt = numel(T);
r  = zeros(1, Nt);
r(1) = r0;
for k = 1:Nt-1
    u = rand;
    if u <= P(r(k),1)
        r(k+1) = 1;
    else
        r(k+1) = 2;
    end
end
end

function Xd = get_delayed_followers(td, T, X, Xhist)
if td <= 0
    Xd = Xhist; return;
end
[i0,i1,w] = uniform_delay_indices(td,T);
Xd = (1-w)*X(:,:,i0) + w*X(:,:,i1);
end

function Sd = get_delayed_leader(td, T, S, Shist)
if td <= 0
    Sd = Shist; return;
end
[i0,i1,w] = uniform_delay_indices(td,T);
Sd = (1-w)*S(:,i0) + w*S(:,i1);
end

function [i0,i1,w] = uniform_delay_indices(td,T)
dt = T(2)-T(1);
pos = (td-T(1))/dt + 1;
nearest = round(pos);
if abs(pos-nearest) < 1e-10
    pos = nearest;
end
pos = max(1,min(pos,numel(T)));
i0 = floor(pos);
i1 = ceil(pos);
w = pos-i0;
end
