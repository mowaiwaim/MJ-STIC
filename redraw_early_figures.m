function redraw_early_figures(outputDir, sourceRoot)
% Redraw the four network simulation plots (Figures 2--5).
% With no arguments, use the accompanying early_figure_data.mat and write
% PDFs alongside this file. To recreate the saved samples, also specify the
% original simulation source root as the second argument.
scriptDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outputDir), outputDir = scriptDir; end
if ~isfolder(outputDir), mkdir(outputDir); end
if nargin < 2
    saved = load(fullfile(scriptDir,'early_figure_data.mat'),'data');
    data = saved.data;
else
    data = sampleOriginalFigures(sourceRoot);
    save(fullfile(outputDir,'early_figure_data.mat'),'data');
end
drawNetworkFigures(data,outputDir);
end

function data = sampleOriginalFigures(sourceRoot)
% Original simulation entry points handle initialization, switching and noise.
% Run from a temporary folder so their CSV logs do not overwrite user data.
oldDir = pwd;
oldPath = path;
oldVisible = get(groot,'defaultFigureVisible');
restore = onCleanup(@() restoreSession(oldDir,oldPath,oldVisible)); %#ok<NASGU>
stageDir = tempname;
mkdir(stageDir);
cd(stageDir);
set(groot,'defaultFigureVisible','off');
addpath(fullfile(sourceRoot,'examples','passivity'));
main_sdde_mj_hybrid_triggerlog_theta_compare;
copyfile(fullfile(stageDir,'trigger_episode_log.csv'), ...
    fullfile(sourceRoot,'outputs','example1_trigger_episode_log.csv'));
copyfile(fullfile(stageDir,'actual_impulse_log.csv'), ...
    fullfile(sourceRoot,'outputs','example1_actual_impulse_log.csv'));
ax = figureAxes(1);
[data.passivity.time,data.passivity.norms] = readCurves(ax,3001,15);
data.passivity.mode2 = readIntervals(ax);
markers = findall(ax,'Type','line','Marker','.');
data.passivity.impulses = markers.XData;

addpath(fullfile(sourceRoot,'examples','synchronization'));
main_sdde_sync_shaded_globalrng_triggerlog_theta_compare;
copyfile(fullfile(stageDir,'trigger_episode_log_sync.csv'), ...
    fullfile(sourceRoot,'outputs','example2_trigger_episode_log.csv'));
copyfile(fullfile(stageDir,'actual_impulse_log_sync.csv'), ...
    fullfile(sourceRoot,'outputs','example2_actual_impulse_log.csv'));
ax = figureAxes(2);
[data.sync.time,data.sync.errors] = readCurves(ax,3001,5);
data.sync.mode2 = readIntervals(ax);
ax = figureAxes(3);
signal = findobj(ax,'Type','stair');
data.sync.mode = signal.YData;
assert(isequal(signal.XData,data.sync.time));
markers = findobj(ax,'Type','line','Marker','o');
data.sync.impulses = markers.XData;

addpath(fullfile(sourceRoot,'examples','uncontrolled'));
main_sdde_mj_nocontrol_shaded_globalrng;
ax = figureAxes(2);
[data.uncontrolled.time,data.uncontrolled.disagreement] = readCurves(ax,10001,2);
close all;
end

function ax = figureAxes(number)
f = findall(groot,'Type','figure','Number',number);
ax = findobj(f,'Type','axes');
assert(isscalar(ax));
end

function [time, curves] = readCurves(ax,nSamples,nCurves)
lines = flipud(findobj(ax,'Type','line'));
lines = lines(arrayfun(@(h) numel(h.XData)==nSamples,lines));
assert(numel(lines)==nCurves);
time = lines(1).XData;
curves = zeros(nCurves,nSamples);
for k = 1:nCurves
    assert(isequal(lines(k).XData,time));
    curves(k,:) = lines(k).YData;
end
end

function intervals = readIntervals(ax)
patches = findall(ax,'Type','patch');
intervals = zeros(numel(patches),2);
for k = 1:numel(patches)
    intervals(k,:) = [min(patches(k).XData) max(patches(k).XData)];
end
intervals = sortrows(intervals);
end

function restoreSession(oldDir,oldPath,oldVisible)
cd(oldDir);
path(oldPath);
set(groot,'defaultFigureVisible',oldVisible);
end

function drawNetworkFigures(data,outputDir)
blue = [0.10 0.36 0.69];
red = [0.76 0.16 0.15];
colors = [0.08 0.37 0.74; 0.81 0.46 0.03; 0.79 0.19 0.15; ...
          0.46 0.25 0.68; 0.11 0.47 0.30];

[f,ax] = plotFrame(8.2);
xlim(ax,[0 3]); ylim(ax,[0 10]);
xticks(ax,0:0.5:3); yticks(ax,0:2:10);
drawIntervals(ax,data.passivity.mode2);
styles = {'-','--',':'};
symbols = {'e','\beta','\alpha'};
curves = gobjects(15,1);
labels = cell(15,1);
for group = 1:3
    for node = 1:5
        k = 5*(group-1)+node;
        curves(k) = plot(ax,data.passivity.time,data.passivity.norms(k,:), ...
            'Color',colors(node,:),'LineStyle',styles{group},'LineWidth',1.0);
        labels{k} = sprintf('$\\|%s_%d(t)\\|_2$',symbols{group},node);
    end
end
drawImpulses(ax,data.passivity.impulses,0);
xlabel(ax,'Time'); ylabel(ax,'Node-wise norm');
order = reshape([1:5; 6:10; 11:15],1,[]);
key = plotKey(ax,curves(order),labels(order),5);
key.FontSize = 8.2;
savePlot(f,outputDir,'figure2.pdf');

[f,ax] = plotFrame(6.6);
xlim(ax,[0 3]); ylim(ax,[0.8 2.4]);
xticks(ax,0:0.5:3); yticks(ax,[1 2]);
signal = stairs(ax,data.sync.time,data.sync.mode,'Color',blue,'LineWidth',1.2);
marks = drawImpulses(ax,data.sync.impulses,0.848);
xlabel(ax,'Time');
ylabel(ax,'Mode');
plotKey(ax,[signal marks],{'$r(t)$','Actual impulses'},1);
savePlot(f,outputDir,'figure3.pdf');

[f,ax] = plotFrame(7.6);
xlim(ax,[0 3]); ylim(ax,[0 0.7]);
xticks(ax,0:0.5:3); yticks(ax,0:0.2:0.6);
drawIntervals(ax,data.sync.mode2);
curves = gobjects(5,1);
labels = cell(5,1);
for node = 1:5
    curves(node) = plot(ax,data.sync.time,data.sync.errors(node,:), ...
        'Color',colors(node,:),'LineWidth',1.05);
    labels{node} = sprintf('$\\|e_%d(t)\\|_2$',node);
end
drawImpulses(ax,data.sync.impulses,0);
xlabel(ax,'Time');
ylabel(ax,'Synchronization error norm');
plotKey(ax,curves,labels,1);
savePlot(f,outputDir,'figure4.pdf');

[f,ax] = plotFrame(7.6);
xlim(ax,[0 10]); ylim(ax,[0 1]);
xticks(ax,0:2:10); yticks(ax,0:0.2:1);
leader = plot(ax,data.uncontrolled.time,data.uncontrolled.disagreement(1,:), ...
    'Color',blue,'LineWidth',1.05);
followers = plot(ax,data.uncontrolled.time,data.uncontrolled.disagreement(2,:), ...
    '--','Color',red,'LineWidth',1.05);
xlabel(ax,'Time');
ylabel(ax,'Disagreement');
plotKey(ax,[leader followers],{'$D_L(t)$','$D_F(t)$'},1);
savePlot(f,outputDir,'figure5.pdf');
end

function [f,ax] = plotFrame(height)
f = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 15.8 height],'Renderer','painters');
ax = axes(f,'Units','normalized','Position',[0.105 0.173 0.865 0.774], ...
    'FontName','Arial','FontUnits','points','FontSize',10,'Box','on', ...
    'TickDir','out','LineWidth',0.65,'Layer','top','XColor',[.14 .14 .14], ...
    'YColor',[.14 .14 .14],'TickLength',[.009 .009]);
hold(ax,'on');
end

function drawIntervals(ax,intervals)
limits = ylim(ax);
for k = 1:size(intervals,1)
    bounds = intervals(k,:);
    patch(ax,bounds([1 2 2 1]),limits([1 1 2 2]),[.94 .94 .94], ...
        'EdgeColor','none','HandleVisibility','off');
end
end

function h = drawImpulses(ax,times,level)
h = plot(ax,times,level*ones(size(times)),'ko','LineStyle','none', ...
    'MarkerSize',3.2,'MarkerFaceColor','k','Clipping','off');
end

function key = plotKey(ax,handles,labels,columnCount)
key = legend(ax,handles,labels,'Interpreter','latex','FontSize',8.5, ...
    'NumColumns',columnCount,'Location','northeast','Box','on', ...
    'EdgeColor',[.75 .75 .75],'Color','white');
key.ItemTokenSize = [16 7];
end

function savePlot(f,outputDir,name)
exportgraphics(f,fullfile(outputDir,name),'ContentType','vector', ...
    'BackgroundColor','white');
close(f);
end
