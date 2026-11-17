function plot_comparison_figure(traj, filename)
% Draw the instantaneous and accumulated errors for the same six schedules.
keys = {'STIC','Zhan','Shang','Liu','PeriodicAll','LuNodeBudget'};
labels = {'MJ-STIC','Zhan et al. (2025)','Shang and Zhang (2025)', ...
          'Liu et al. (2026)','Periodic (all nodes)', ...
          'Lu et al. (2012)'};
% Okabe-Ito palette, excluding yellow against the white background.
colors = [0 114 178; 213 94 0; 0 158 115; 204 121 167; ...
          230 159 0; 86 180 233] / 255;
styles = {'-','--',':','-.','-.','--'};
panelLabels = {'(a)','(b)'};
t = traj.Time;
errors = traj{:,keys};
integrated = cumtrapz(t,errors);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 18.8 10.6], ...
             'Color','white','Renderer','painters');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
tiles = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');
axesList = gobjects(1,2);
for panel = 1:2
    ax = nexttile(tiles);
    axesList(panel) = ax;
    hold(ax,'on');
    if panel == 1
        values = errors;
        set(ax,'YScale','log','YLim',[1e-10 1], 'YTick',10.^(-10:2:0));
        ylabel(ax,'Mean-square error');
    else
        values = integrated;
        ylim(ax,[0 .3]);
        yticks(ax,0:.1:.3);
        ylabel(ax,'Integrated mean-square error');
    end
    handles = gobjects(1,numel(keys));
    for k = 1:numel(keys)
        handles(k) = plot(ax,t,values(:,k),'LineStyle',styles{k}, ...
            'LineWidth',1.1 + .5*(k==1),'Color',colors(k,:), ...
            'DisplayName',labels{k});
    end
    set(ax,'FontName','Arial','FontSize',9,'TickDir','out','Box','off', ...
        'XLim',[0 t(end)],'XTick',0:2:t(end),'TickLength',[.013 .015]);
    xlabel(ax,'Time');
    title(ax,panelLabels{panel},'FontSize',10,'FontWeight','normal');
end
lg = legend(axesList(2),handles,labels,'Orientation','horizontal', ...
    'NumColumns',3,'Box','off','FontName','Arial','FontSize',8);
lg.ItemTokenSize = [23 9];
lg.Layout.Tile = 'south';

% Show the first impulses at a linear scale within the unused lower left
% of panel (a). The main curves and their full time range remain visible.
drawnow;
pos = axesList(1).InnerPosition;
inset = axes('Parent',fig,'Units','normalized', ...
    'Position',[pos(1)+.15*pos(3),pos(2)+.12*pos(4),.47*pos(3),.27*pos(4)], ...
    'FontName','Arial','FontSize',7.1,'Color','white','Box','on','TickDir','out');
hold(inset,'on');
early = t <= .6;
for k = numel(keys):-1:1
    plot(inset,t(early),errors(early,k),'LineWidth',.85+.4*(k==1), ...
        'Color',colors(k,:),'LineStyle',styles{k});
end
set(inset,'XLim',[0 .6],'XTick',[0 .3 .6],'YLim',[0 .5],'YTick',[0 .25 .5], ...
    'TickLength',[.035 .04]);
exportgraphics(fig,filename,'ContentType','vector','BackgroundColor','white');
end
