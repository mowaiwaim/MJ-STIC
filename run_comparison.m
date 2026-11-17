function run_comparison()
% Paired comparison with the largest-error pinning law of Lu et al. (2012),
% IEEE TNNLS 23(2), 285-292, Eqs. (6)-(7), doi:10.1109/TNNLS.2011.2179312.
% Run from this folder. The same CTMC and Brownian paths drive all methods.
p = bw_parameters();
outdir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~isfolder(outdir), mkdir(outdir); end
nodeCount = numel(p.x0);

% Tune only pulse counts on an independent pilot set; test errors are never
% consulted. Candidate periods are specified in physical time (seconds).
pilot = p;
pilot.seed = 512091;
pilot.samples = 256;
trainingPaths = bw_paths(pilot);
training = simulate_bounded_wait(pilot,trainingPaths,'stic',NaN);
pilotEvents = mean(training.impulses);
testPeriods = (50:1000)'/1000;
dwell = restricted_episodes(trainingPaths,pilot.T);
eventsAtPeriod = arrayfun(@(d) sum(floor(dwell./d+1e-12))/pilot.samples, testPeriods);
targets = [pilotEvents,nodeCount*pilotEvents];
chosen = zeros(2,1);
for k=1:2
    assert(min(eventsAtPeriod)<targets(k) && max(eventsAtPeriod)>targets(k));
    [~,idx] = min(abs(eventsAtPeriod-targets(k)));
    chosen(k) = testPeriods(idx);
end
pilotTable = table(["Event budget";"Node correction budget"],targets',chosen, ...
    arrayfun(@(d) sum(floor(dwell./d+1e-12))/pilot.samples,chosen), ...
    'VariableNames',{'Budget','TargetImpulses','Period','PilotImpulses'});
writetable(pilotTable,fullfile(outdir,'comparison_pilot.csv'));

methods = ["MJ-STIC";"Periodic all-node"; ...
           "Lu pinning, node budget"];
codes = {'stic','periodic','lu'};
intervals = [NaN;chosen(1);chosen(2)];
responses = cell(numel(codes),1);
testPaths = bw_paths(p);
summary = table(methods,intervals,'VariableNames',{'Method','Period'});
pathRows = cell(numel(codes),1);
for k=1:numel(codes)
    sim = simulate_bounded_wait(p,testPaths,codes{k},intervals(k));
    responses{k} = sim;
    corrections = sim.impulses;
    if ~strcmp(codes{k},'lu'), corrections = nodeCount*corrections; end
    % One snapshot means one full-vector error acquisition in the constrained
    % mode. The self-trigger computes Q only on opening a reference. An
    % executed pulse also opens the next reference and reuses that sample.
    % Fixed clocks need the full error only when executing a pulse: Lu's rule
    % ranks all nodes, the all-node control uses each node's error.
    snapshots = sim.impulses;
    if strcmp(codes{k},'stic'), snapshots = sim.references; end
    summary.TerminalMSE(k) = mean(sim.terminal);
    summary.IntegratedMSE(k) = mean(sim.integral);
    summary.Impulses(k) = mean(sim.impulses);
    summary.NodeCorrections(k) = mean(corrections);
    summary.FullErrorSnapshots(k) = mean(snapshots);
    summary.CancelledTimers(k) = mean(sim.cancellations);
    summary.TimerStarts(k) = mean(sim.references);
    pathRows{k} = table((1:p.samples)',repmat(methods(k),p.samples,1),sim.terminal, ...
        sim.integral,sim.impulses,corrections,snapshots,sim.cancellations,sim.references, ...
        'VariableNames',{'Trial','Method','TerminalSquaredError','IntegralSquaredError', ...
        'Impulses','NodeCorrections','FullErrorSnapshots','CancelledTimers','TimerStarts'});
end
writetable(summary,fullfile(outdir,'comparison_summary.csv'));
writetable(vertcat(pathRows{:}),fullfile(outdir,'comparison_paths.csv'));

% Paired percentile bootstrap confidence intervals across independent
% realizations. Positive differences mean greater error than our method.
bootStream = RandStream('mt19937ar','Seed',39147);
bootIdx = randi(bootStream,p.samples,p.samples,5000);
contrasts = table();
for k=2:numel(codes)
    for item = {'terminal','integral'}
        key = item{1};
        delta = responses{k}.(key) - responses{1}.(key);
        ordered = sort(mean(delta(bootIdx),1));
        interval = ordered([125,4875]);
        contrasts = [contrasts; table(methods(k),string(key),mean(delta), ...
            interval(1),interval(2),'VariableNames', ...
            {'MethodMinusSTIC','Metric','MeanDiff','Lower95','Upper95'})]; %#ok<AGROW>
    end
end
writetable(contrasts,fullfile(outdir,'comparison_paired.csv'));
traj = table(responses{1}.time,'VariableNames',{'Time'});
curveNames = {'STIC','PeriodicAll','LuNodeBudget'};
for k=1:numel(codes), traj.(curveNames{k})=mean(responses{k}.energy,2); end
writetable(traj,fullfile(outdir,'comparison_trajectories.csv'));
add_recent_comparisons(p,outdir,bootIdx);
disp(readtable(fullfile(outdir,'comparison_summary.csv')));
end

function dwell = restricted_episodes(paths,T)
% The exact CTMC starts in mode 2 and alternates modes. A fixed-period
% timer starts on each mode-2 entry and stops on exit or at the horizon.
parts = cellfun(@(times) diff([0,times,T]),paths.switches,'UniformOutput',false);
parts = cellfun(@(v) v(1:2:end)',parts,'UniformOutput',false);
dwell = vertcat(parts{:});
end
