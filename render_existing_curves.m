function render_existing_curves
folder=fileparts(mfilename('fullpath'));
p=bw_parameters(); paths=bw_paths(p);
curve=readtable(fullfile(folder,'bw_mean_curves.csv'));
M=[curve.stic,curve.zero_floor,curve.noImpulses];
sample=readtable(fullfile(folder,'bw_first_path.csv'));
proposed.time=sample.time;
proposed.energy=sample.squaredError;
proposed.firstEvents=readmatrix(fullfile(folder,'bw_first_events.csv'));
assert(isequal(curve.time,sample.time));
fig=figure('Visible','off','Color','white','Units','inches','Position',[1 1 8.0 3.20]);
closeFig=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);
hold(ax,'on'); box(ax,'on');
% Plot exact pre/post impulse values between grid samples.
tMax=p.T; if ~isempty(paths.switches{1}), tMax=paths.switches{1}(1); end
idx=find(proposed.time<=tMax+1e-10);
t=proposed.time(idx); energy=proposed.energy(idx,1);
events=proposed.firstEvents(proposed.firstEvents(:,1)<=tMax,:);
tAll=[t;events(:,1);events(:,1)];
eAll=[energy;events(:,2);events(:,3)];
[~,ord]=sortrows([tAll,[3*ones(numel(t),1);ones(size(events,1),1);2*ones(size(events,1),1)]],[1 2]);
assert(isempty(paths.switches{1}) || paths.switches{1}(1)>=tMax, ...
    'Panel (a) uses the initial uninterrupted constrained-mode episode.');
palette=[0 114 178;213 94 0;0 158 115]/255;
plot(ax,tAll(ord),eAll(ord),'Color',palette(1,:),'LineWidth',1.4);
plot(ax,events(:,1),events(:,2),'o','Color',palette(2,:),'MarkerFaceColor','white','MarkerSize',4.8,'LineWidth',1.1);
plot(ax,events(:,1),events(:,3),'o','Color',palette(2,:),'MarkerFaceColor',palette(2,:),'MarkerSize',4.8);
xlim(ax,[0 tMax]); ylim(ax,[0 0.48]);
xlabel(ax,'Time'); ylabel(ax,'Squared error');
title(ax,'(a)','FontWeight','normal');
ax=nexttile(layout); hold(ax,'on'); box(ax,'on'); 
set(ax,'YScale','log'); styles={'-',':','-.'};
h=gobjects(4,1);
for k=1:3
    h(k)=plot(ax,proposed.time,M(:,k),'LineWidth',1.4,'Color',palette(k,:),'LineStyle',styles{k});
end
h(4)=plot(ax,proposed.time,curve.certificate,'--','Color',[0.2 0.2 0.2],'LineWidth',1.2);
xlabel(ax,'Time'); ylabel(ax,'Mean-square error');
title(ax,'(b)','FontWeight','normal'); xlim(ax,[0 p.T]);
legend(ax,h,{'MJ-STIC','No maximum waiting time','No control in mode {\bf 2}','Theoretical error bound'}, ...
    'Location','southwest','Box','off','FontSize',8);
set(findall(fig,'Type','axes'),'FontName','Arial','FontSize',9,'TickDir','out','Box','off');
set(findall(fig,'Type','legend'),'FontName','Arial');
exportgraphics(fig,fullfile(folder,'figures','figure6.pdf'),'ContentType','vector','BackgroundColor','white');
end
