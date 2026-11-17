function main_sdde_mj_hybrid_triggerlog_theta_compare
% Markov jump SDDE leader-follower network with
% mode 1: implicit continuous leader-pinning PD control only
% mode 2: self-triggered impulsive control only
%
% This version uses the revised my7paper dynamics and adds trigger logs:
%   1) episode-level log for each reference instant theta
%   2) actual-impulse log
%   3) threshold Ptheta_eps: if e(theta)^T M e(theta) <= Ptheta_eps,
%      then hat_t(theta) is treated as Inf on the simulation grid.
%
% Output files written to the current MATLAB folder when you run:
%   trigger_episode_log.csv
%   actual_impulse_log.csv
% Figure 1 is a node-wise diagnostic; Figure 2 is the three-curve node-RMS
% publication figure exported by export_my8paper_figures.m; Figure 3 is a
% passivity-mode/impulse diagnostic and is not used as the paper mode plot.

clc; close all;
clear new_sddefun_mj;
rng(2,'twister');

%% ===================== 1) basic settings =====================
t0     = 0;
tfinal = 3;
dt     = 1e-3;
T      = t0:dt:tfinal;
Nt     = numel(T);

n = 3;
N = 5;

% Use the same constant history as the synchronization/open-loop comparison.
Xhist = [ 0.02 -0.03  0.04 -0.05  0.06;
          0.01 -0.01  0.01 -0.02  0.02;
          0.05 -0.07  0.10 -0.13  0.16];
Shist = [0.3;-0.2;0.4];
X = zeros(n, N, Nt);
S = zeros(n, Nt);
X(:,:,1) = Xhist;
S(:,1)   = Shist;

mu2        = -0.2;
lambda2    = 8.0;
betahat2   = 0.05;
alpha2_st  = 0.3;
Ptheta_eps = 1e-10;   % user can change; if P(theta)<=eps => hat_t = Inf

% Same Markov path as the synchronization example (same T, dt, r0, seed).
Gamma = [-2 2; 8 -8];
r0    = 2;  % start in the impulsive mode so theta=0 is exercised
r     = simulate_markov_chain_grid(Gamma, T, r0);
fprintf('Markov path: %d switches, mode-2 fraction = %.4f\n', ...
    nnz(diff(r)),mean(r==2));

Mmat = get_trigger_matrix();

impulse_flag = false(1, Nt);
next_imp_idx = inf;
theta_idx    = NaN;
theta_active = false;
current_episode_id = NaN;
pending_suspend_episode = NaN;

%% ===== logs =====
log_theta          = [];
log_Ptheta         = [];
log_hat_t          = [];
log_zeta           = [];
log_compare        = strings(0,1);
log_action         = strings(0,1);
log_source         = strings(0,1);
log_theta_next     = [];
log_next_rule      = strings(0,1);
log_mode_at_theta  = [];
log_eps            = [];

imp_theta          = [];
imp_time           = [];
imp_Ptheta         = [];
imp_zeta           = [];
imp_source         = strings(0,1);
imp_episode_id     = [];

% if initially in impulsive mode, theta = 0
if r(1) == 2
    theta_active = true;
    theta_idx    = 1;
    e0 = reshape(X(:,:,1) - repmat(S(:,1), 1, N), [], 1);
    [next_imp_idx, zeta_idx, Ptheta, cmpStr, actionStr, nextRule] = ...
        compute_next_impulse_info(theta_idx, e0, T, r, Mmat, lambda2, betahat2, alpha2_st, Ptheta_eps);

    current_episode_id = numel(log_theta) + 1;
    [log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, log_source, ...
        log_theta_next, log_next_rule, log_mode_at_theta, log_eps] = ...
        append_episode_log(log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, ...
        log_source, log_theta_next, log_next_rule, log_mode_at_theta, log_eps, ...
        T(theta_idx), Ptheta, idx_to_time(T,next_imp_idx), idx_to_time(T,zeta_idx), cmpStr, actionStr, ...
        "initial_in_I", NaN, nextRule, r(theta_idx), Ptheta_eps);
end

%% ===================== 2) Euler-Maruyama time stepping =====================
for k = 1:Nt-1
    tk = T(k);

    hk  = new_delay(tk, []);
    tdk = tk - hk;

    Xd = get_delayed_followers(tdk, T, X, Xhist);
    Sd = get_delayed_leader(tdk, T, S, Shist);

    Alpha = alpha_input_local(tk, N);

    [driftX, driftS, sigmaX, sigmaS] = new_sddefun_mj(tk, X(:,:,k), Xd, S(:,k), Sd, r(k), Alpha);

    dW = sqrt(dt) * randn;

    X(:,:,k+1) = X(:,:,k) + driftX * dt + sigmaX * dW;
    S(:,k+1)   = S(:,k)   + driftS * dt + sigmaS * dW;

    % entering impulsive mode: theta := entry time into I={2}
    if r(k) ~= 2 && r(k+1) == 2
        theta_active = true;
        theta_idx    = k+1;
        e_theta = reshape(X(:,:,k+1) - repmat(S(:,k+1), 1, N), [], 1);
        [next_imp_idx, zeta_idx, Ptheta, cmpStr, actionStr, nextRule] = ...
            compute_next_impulse_info(theta_idx, e_theta, T, r, Mmat, lambda2, betahat2, alpha2_st, Ptheta_eps);

        current_episode_id = numel(log_theta) + 1;
        [log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, log_source, ...
            log_theta_next, log_next_rule, log_mode_at_theta, log_eps] = ...
            append_episode_log(log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, ...
            log_source, log_theta_next, log_next_rule, log_mode_at_theta, log_eps, ...
            T(theta_idx), Ptheta, idx_to_time(T,next_imp_idx), idx_to_time(T,zeta_idx), cmpStr, actionStr, ...
            "enter_I", NaN, nextRule, r(theta_idx), Ptheta_eps);

        if ~isnan(pending_suspend_episode)
            log_theta_next(pending_suspend_episode) = T(theta_idx);
            log_next_rule(pending_suspend_episode) = "theta := next entry into I";
            pending_suspend_episode = NaN;
        end
    end

    % leaving impulsive mode: suspend the trigger mechanism
    if r(k) == 2 && r(k+1) ~= 2
        theta_active = false;
        theta_idx    = NaN;
        next_imp_idx = inf;

        if ~isnan(current_episode_id)
            if isinf(log_hat_t(current_episode_id)) || log_hat_t(current_episode_id) >= T(k+1)
                pending_suspend_episode = current_episode_id;
            end
        end
    end

    % actual impulsive instant: t_{k+1} = hat_t(theta) < zeta(theta)
    if theta_active && (k+1 == next_imp_idx) && r(k+1) == 2
        X(:,:,k+1) = X(:,:,k+1) + mu2 * (X(:,:,k+1) - repmat(S(:,k+1), 1, N));
        impulse_flag(k+1) = true;

        % write impulse log for the episode that produced it
        imp_theta(end+1,1)      = log_theta(current_episode_id); %#ok<AGROW>
        imp_time(end+1,1)       = T(k+1); %#ok<AGROW>
        imp_Ptheta(end+1,1)     = log_Ptheta(current_episode_id); %#ok<AGROW>
        imp_zeta(end+1,1)       = log_zeta(current_episode_id); %#ok<AGROW>
        imp_source(end+1,1)     = log_source(current_episode_id); %#ok<AGROW>
        imp_episode_id(end+1,1) = current_episode_id; %#ok<AGROW>

        % theta := t_{k+1} after the impulse
        theta_idx = k+1;
        e_theta = reshape(X(:,:,k+1) - repmat(S(:,k+1), 1, N), [], 1);
        [next_imp_idx, zeta_idx, Ptheta, cmpStr, actionStr, nextRule] = ...
            compute_next_impulse_info(theta_idx, e_theta, T, r, Mmat, lambda2, betahat2, alpha2_st, Ptheta_eps);

        current_episode_id = numel(log_theta) + 1;
        [log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, log_source, ...
            log_theta_next, log_next_rule, log_mode_at_theta, log_eps] = ...
            append_episode_log(log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, ...
            log_source, log_theta_next, log_next_rule, log_mode_at_theta, log_eps, ...
            T(theta_idx), Ptheta, idx_to_time(T,next_imp_idx), idx_to_time(T,zeta_idx), cmpStr, actionStr, ...
            "post_impulse", NaN, nextRule, r(theta_idx), Ptheta_eps);
    end
end

%% ===================== 3) post-processing =====================
E     = zeros(n, N, Nt);
Beta  = zeros(4, N, Nt);
AlphaAll = zeros(n, N, Nt);

E1 = [1.0 0   0;
      0   0.9 0;
      1.1 0   0.8;
      0   1.2 0.7];

E2 = [1.2 0   0;
      0   0.9 1.1;
      0   0   1.0;
      1.2 0   0];

for k = 1:Nt
    AlphaAll(:,:,k) = alpha_input_local(T(k), N);
    E(:,:,k) = X(:,:,k) - repmat(S(:,k), 1, N);
    for i = 1:N
        Beta(:,i,k) = E1 * E(:,i,k) + E2 * AlphaAll(:,i,k);
    end
end

nodeErrorNorm = squeeze(sqrt(sum(E.^2,1)));
nodeOutputNorm = squeeze(sqrt(sum(Beta.^2,1)));
nodeInputNorm = squeeze(sqrt(sum(AlphaAll.^2,1)));

% Paper quantity: RMS over the five node-wise Euclidean norms.
errorRms = sqrt(mean(nodeErrorNorm.^2,1));
outputRms = sqrt(mean(nodeOutputNorm.^2,1));
inputRms = sqrt(mean(nodeInputNorm.^2,1));

mode2_intervals = get_mode2_intervals(T, r);
impulse_times   = T(impulse_flag);

%% ===================== 4) plots =====================
colors = [0.12 0.56 1.00;
          1.00 0.84 0.00;
          0.90 0.25 0.25;
          0.60 0.20 0.98;
          0.18 0.55 0.34];

% Figure 1 is a node-wise diagnostic and is not the paper figure.
figure(1); clf; hold on; grid on;
set(gcf,'Name','Diagnostic: passivity node-wise norms','NumberTitle','off');
for i = 1:N
    plot(T,nodeErrorNorm(i,:),'-','Color',colors(i,:),'LineWidth',1.2);
end
for i = 1:N
    plot(T,nodeOutputNorm(i,:),'--','Color',colors(i,:),'LineWidth',1.2);
end
for i = 1:N
    plot(T,nodeInputNorm(i,:),':','Color',colors(i,:),'LineWidth',1.0);
end
xlim([T(1) T(end)]); ylim([0 10]);
xlabel('Time (s)'); ylabel('Node-wise norm');
leg = cell(1,3*N);
for i = 1:N
    leg{i} = sprintf('$\\|e_%d(t)\\|_2$',i);
    leg{N+i} = sprintf('$\\|\\beta_%d(t)\\|_2$',i);
    leg{2*N+i} = sprintf('$\\|\\alpha_%d(t)\\|_2$',i);
end
legend(leg,'Interpreter','latex','FontSize',9,'Location','northeast', ...
    'NumColumns',3);
add_mode2_shading_and_xaxis_impulses(gca,mode2_intervals,impulse_times);
ax = gca;
ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.0;
hold off;

% Figure 2 is the unique passivity figure used by the paper.
figure(2); clf; hold on; grid on;
set(gcf,'Name','Paper figure: passivity node-RMS norms','NumberTitle','off');
plot(T,errorRms,'-','Color',[0.05 0.35 0.75],'LineWidth',1.2);
plot(T,outputRms,'--','Color',[0.80 0.15 0.15],'LineWidth',1.2);
plot(T,inputRms,':','Color',[0.15 0.55 0.25],'LineWidth',1.2);
xlim([T(1) T(end)]);
ymax = max([errorRms,outputRms,inputRms]);
ylim([0 max(1,1.05*ymax)]);
xlabel('Time (s)'); ylabel('Node RMS norm');
legend({'$e_{\mathrm{NRMS}}$', ...
        '$\beta_{\mathrm{NRMS}}$', ...
        '$\alpha_{\mathrm{NRMS}}$'}, ...
    'Interpreter','latex','FontSize',8,'Location','best', ...
    'Orientation','horizontal');
box on;
hold off;

figure(3); clf; hold on; grid on;

hMode = stairs(T, r, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.1);

xlabel('t(s)');
ylim([0.8 2.4]);
yticks([1 2]);

% 把脉冲点画成黑色，并放在x轴附近
if any(impulse_flag)
    yl = ylim(gca);
    y_imp = yl(1) + 0.03*(yl(2)-yl(1));   % 略高于x轴，避免被边框压住
    hImpulse = plot(T(impulse_flag), y_imp*ones(1,nnz(impulse_flag)), ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    hleg = legend([hMode hImpulse], {'r(t)','actual impulses'}, ...
        'FontSize', 12, 'Location', 'northeast');
else
    hleg = legend(hMode, {'r(t)'}, 'FontSize', 12, 'Location', 'northeast');
end

hleg.NumColumns = 1;

% 坐标轴矩形边框全黑
ax = gca;
ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.0;

hold off;
%% ===================== 5) trigger data output =====================
EpisodeID = (1:numel(log_theta)).';
trigger_episode_table = table(EpisodeID, log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, ...
    log_action, log_source, log_theta_next, log_next_rule, log_mode_at_theta, log_eps, ...
    'VariableNames', {'EpisodeID','theta','Ptheta','hat_t','zeta','compare_hat_vs_zeta', ...
    'action','theta_source','theta_next_actual','theta_update_rule','mode_at_theta','Ptheta_eps'});

ImpulseID = (1:numel(imp_time)).';
actual_impulse_table = table(ImpulseID, imp_episode_id, imp_theta, imp_time, imp_Ptheta, imp_zeta, imp_source, ...
    'VariableNames', {'ImpulseID','EpisodeID','theta','t_k','Ptheta','zeta','theta_source'});

fprintf('\n==================== Trigger episode log ====================\n');
disp(trigger_episode_table);

fprintf('\n==================== Actual impulse log ====================\n');
disp(actual_impulse_table);

writetable(trigger_episode_table, 'trigger_episode_log.csv');
writetable(actual_impulse_table, 'actual_impulse_log.csv');

fprintf('\nCSV files written:\n');
fprintf('  %s\n', fullfile(pwd, 'trigger_episode_log.csv'));
fprintf('  %s\n', fullfile(pwd, 'actual_impulse_log.csv'));
fprintf('Total actual impulsive instants: %d\n', nnz(impulse_flag));

end

%% ===================== local helper functions =====================
function [next_imp_idx, zeta_idx, Ptheta, cmpStr, actionStr, nextRule] = ...
    compute_next_impulse_info(theta_idx, e_theta, T, r, Mmat, lambda2, betahat2, alpha2_st, Ptheta_eps)
% Compute the candidate impulsive instant \hat t(theta) on the whole
% remaining simulation horizon first, and then compare it with zeta(theta).
% If Ptheta <= Ptheta_eps, we treat hat_t(theta)=Inf.

next_imp_idx = inf;
zeta_idx = inf;
Ptheta = NaN;
cmpStr = "invalid";
actionStr = "invalid theta";
nextRule = "none";

if isnan(theta_idx) || theta_idx >= numel(T) || r(theta_idx) ~= 2
    return;
end

Ptheta = e_theta.' * Mmat * e_theta;

zeta_rel = find(r(theta_idx+1:end) ~= r(theta_idx), 1, 'first');
if ~isempty(zeta_rel)
    zeta_idx = theta_idx + zeta_rel;
end

if Ptheta <= Ptheta_eps
    next_imp_idx = inf;
    cmpStr = "hat_t = Inf (P small)";
    if isfinite(zeta_idx)
        actionStr = "discard candidate, suspend at zeta";
        nextRule = "theta suspended; restart when re-enter I";
    else
        actionStr = "no trigger on grid horizon";
        nextRule = "none on current horizon";
    end
    return;
end

% IMPORTANT: candidate time is computed on the whole remaining horizon,
% not truncated to the current mode-2 interval.
for j = theta_idx+1 : numel(T)
    tau = T(j) - T(theta_idx);
    xi  = Ptheta * (exp(lambda2 * tau) - 1) - betahat2 * Ptheta - exp(-alpha2_st * tau);
    if xi >= 0
        next_imp_idx = j;
        break;
    end
end

if isfinite(next_imp_idx)
    if isinf(zeta_idx) || next_imp_idx < zeta_idx
        cmpStr = "hat_t < zeta";
        actionStr = "actual impulse at hat_t";
        nextRule = "theta := hat_t after impulse";
    else
        cmpStr = "hat_t >= zeta";
        actionStr = "discard candidate, suspend at zeta";
        nextRule = "theta suspended; restart when re-enter I";
    end
else
    if isfinite(zeta_idx)
        cmpStr = "hat_t = Inf >= zeta";
        actionStr = "discard candidate, suspend at zeta";
        nextRule = "theta suspended; restart when re-enter I";
    else
        cmpStr = "hat_t = Inf on current horizon";
        actionStr = "no trigger on grid horizon";
        nextRule = "none on current horizon";
    end
end
end

function [log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, log_source, ...
    log_theta_next, log_next_rule, log_mode_at_theta, log_eps] = ...
    append_episode_log(log_theta, log_Ptheta, log_hat_t, log_zeta, log_compare, log_action, ...
    log_source, log_theta_next, log_next_rule, log_mode_at_theta, log_eps, ...
    theta, Ptheta, hat_t, zeta, cmpStr, actionStr, sourceStr, thetaNext, nextRule, modeVal, epsVal)

log_theta(end+1,1)         = theta;
log_Ptheta(end+1,1)        = Ptheta;
log_hat_t(end+1,1)         = hat_t;
log_zeta(end+1,1)          = zeta;
log_compare(end+1,1)       = string(cmpStr);
log_action(end+1,1)        = string(actionStr);
log_source(end+1,1)        = string(sourceStr);
log_theta_next(end+1,1)    = thetaNext;
log_next_rule(end+1,1)     = string(nextRule);
log_mode_at_theta(end+1,1) = modeVal;
log_eps(end+1,1)           = epsVal;
end

function tval = idx_to_time(T, idx)
if isinf(idx) || isnan(idx)
    tval = inf;
else
    tval = T(idx);
end
end

function add_mode2_shading_and_xaxis_impulses(ax, mode2_intervals, impulse_times)
yl = ylim(ax);
for k = 1:size(mode2_intervals,1)
    x1 = mode2_intervals(k,1);
    x2 = mode2_intervals(k,2);
    p = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(p, 'bottom');
    xline(ax, x1, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, x2, ':', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
if ~isempty(impulse_times)
    plot(ax, impulse_times, yl(1)*ones(size(impulse_times)), 'k.', ...
        'MarkerSize', 15, 'HandleVisibility', 'off', 'Clipping', 'off');
end
ylim(ax, yl);
end

function intervals = get_mode2_intervals(T, r)
mask = (r == 2);
starts = find(diff([0 mask]) == 1);
stops  = find(diff([mask 0]) == -1);
intervals = zeros(numel(starts), 2);
for k = 1:numel(starts)
    i1 = starts(k);
    i2 = stops(k);
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
Xd = zeros(size(Xhist));
for i = 1:size(X,2)
    for p = 1:size(X,1)
        Xd(p,i) = interp1(T, squeeze(X(p,i,:)), td, 'linear', 'extrap');
    end
end
end

function Sd = get_delayed_leader(td, T, S, Shist)
if td <= 0
    Sd = Shist; return;
end
Sd = zeros(size(Shist));
for p = 1:size(S,1)
    Sd(p) = interp1(T, S(p,:), td, 'linear', 'extrap');
end
end

function Mmat = get_trigger_matrix()
n = 3; N = 5;
[Psi,~,W,bD,~] = hybrid_params();
Mmat = eye(n*N);
for q = 1:3
    Mmat = Mmat - bD(q) * kron(W{q}, Psi{q});
end
end

function [Psi,H,W,bD,bP] = hybrid_params()
Psi = cell(1,3);
Psi{1} = diag([0.2 0.3 0.6]);
Psi{2} = diag([0.8 0.2 0.4]);
Psi{3} = diag([0.3 0.4 0.5]);

H = cell(1,3);
H{1} = [-0.3 0.1 0.2 0   0; 0.1 -0.2 0 0.1 0; 0.2 0 -0.4 0 0.2; 0 0.1 0 -0.3 0.2; 0 0 0.2 0.2 -0.4];
H{2} = [-0.4 0.2 0.2 0   0; 0.2 -0.4 0 0.2 0; 0.2 0 -0.3 0 0.1; 0 0.2 0 -0.3 0.1; 0 0 0.1 0.1 -0.2];
H{3} = [-0.5 0.2 0.3 0   0; 0.2 -0.4 0 0.2 0; 0.3 0 -0.5 0 0.2; 0 0.2 0 -0.3 0.1; 0 0 0.2 0.1 -0.3];

W = cell(1,3);
W{1} = [-0.2 0.1 0.1 0   0; 0.1 -0.2 0 0.1 0; 0.1 0 -0.3 0 0.2; 0 0.1 0 -0.3 0.2; 0 0 0.2 0.2 -0.4];
W{2} = [-0.3 0.1 0.2 0   0; 0.1 -0.3 0 0.2 0; 0.2 0 -0.3 0 0.1; 0 0.2 0 -0.3 0.1; 0 0 0.1 0.1 -0.2];
W{3} = [-0.3 0.2 0.1 0   0; 0.2 -0.3 0 0.1 0; 0.1 0 -0.3 0 0.2; 0 0.1 0 -0.2 0.1; 0 0 0.2 0.1 -0.3];

bD = [0.3 0.4 0.5];
bP = [15.2387 17.5570 15.8809];
end

function Alpha = alpha_input_local(t, N)
n = 3;
Alpha = zeros(n, N);
amp = (1:N) / 2;
s   = abs(sin(pi * t));
for i = 1:N
    Alpha(:,i) = amp(i) * s * ones(n,1);
end
end
