function add_recent_comparisons(p, outdir, bootIdx)
% Add the three adapted self-triggered rules to the same held-out comparison.
sourceDir = fileparts(mfilename('fullpath'));
run_zhan_comparison();
run_shang_comparison();
run_liu_static_comparison(sourceDir,outdir);
writetable(verify_liu_static_certificate(p),fullfile(outdir,'liu_certificate.csv'));

summaryFile = fullfile(outdir,'comparison_summary.csv');
pathFile = fullfile(outdir,'comparison_paths.csv');
curveFile = fullfile(outdir,'comparison_trajectories.csv');
pairedFile = fullfile(outdir,'comparison_paired.csv');
summary = readtable(summaryFile,'TextType','string');
paths = readtable(pathFile,'TextType','string');
traj = readtable(curveFile);
contrasts = readtable(pairedFile,'TextType','string');
summary.ScalarErrorSamples = numel(p.x0)*summary.FullErrorSnapshots;
paths.ScalarErrorSamples = numel(p.x0)*paths.FullErrorSnapshots;
reference = paths(paths.Method=="MJ-STIC",:);
names = ["Zhan et al. (2025)";"Shang and Zhang (2025)";"Liu et al. (2026)"];
keys = ["zhan";"shang";"liu"];
curveKeys = ["Zhan";"Shang";"Liu"];
assert(height(reference)==p.samples && isequal(reference.Trial,(1:p.samples)'));
for k=1:numel(keys)
    key = keys(k);
    record = readtable(fullfile(outdir,key+"_paths.csv"));
    curve = readtable(fullfile(outdir,key+"_trajectories.csv"));
    assert(isequal(record.Trial,reference.Trial) && isequal(curve.Time,traj.Time));
    assert(abs(mean(record.TerminalSquaredError)-curve{end,2}) < 1e-12);
    assert(abs(mean(record.IntegralSquaredError)-trapz(curve.Time,curve{:,2})) < 1e-10);
    assert(all(isfinite(record.IntegralSquaredError)) && all(record.NodeCorrections>=record.Impulses));
    n=height(record);
    if key=="liu"
        scalars=record.ScalarErrorSamples;
    else
        scalars=numel(p.x0)*record.FullErrorSnapshots;
    end
    % Timer bookkeeping differs for multistage and asynchronous schedules.
    % The individual method files retain their own cancellation/probe counts.
    unspecified=NaN(n,1);
    rows=table(record.Trial,repmat(names(k),n,1),record.TerminalSquaredError, ...
        record.IntegralSquaredError,record.Impulses,record.NodeCorrections, ...
        record.FullErrorSnapshots,unspecified,unspecified,scalars, ...
        'VariableNames',paths.Properties.VariableNames);
    paths=[paths;rows]; %#ok<AGROW>
    avg=table(names(k),NaN,mean(record.TerminalSquaredError), ...
        mean(record.IntegralSquaredError),mean(record.Impulses), ...
        mean(record.NodeCorrections),mean(record.FullErrorSnapshots),NaN,NaN, ...
        mean(scalars),'VariableNames',summary.Properties.VariableNames);
    summary=[summary;avg]; %#ok<AGROW>
    traj.(curveKeys(k))=curve{:,2};
    metrics=["terminal","integral"];
    fields=["TerminalSquaredError","IntegralSquaredError"];
    for m=1:2
        delta=record.(fields(m))-reference.(fields(m));
        boot=sort(mean(delta(bootIdx),1));
        entry=table(names(k),metrics(m),mean(delta),boot(125),boot(4875), ...
            'VariableNames',contrasts.Properties.VariableNames);
        contrasts=[contrasts;entry]; %#ok<AGROW>
    end
end
writetable(summary,summaryFile);
writetable(paths,pathFile);
writetable(traj,curveFile);
writetable(contrasts,pairedFile);
plot_comparison_figure(traj,fullfile(fileparts(outdir),'figures','figure7.pdf'));
end
