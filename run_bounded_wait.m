function run_bounded_wait()
% MATLAB R2024a/base MATLAB. Run from this folder, or add it to the path.
folder=fileparts(mfilename('fullpath'));
p=bw_parameters();
cert=bw_certificate(p);
paths=bw_paths(p);
names={'stic','zero_floor','none'};
displayNames=["MJ-STIC";"zero_floor";"none"];
runs=cell(1,numel(names));
for k=1:numel(names)
    runs{k}=simulate_bounded_wait(p,paths,names{k},NaN);
end
proposed=runs{1};
assert(all(proposed.maxAge <= p.H+1e-10));
assert(all(abs(proposed.firstEvents(:,3)-(1+p.mu)^2*proposed.firstEvents(:,2))<1e-12));
nt=numel(proposed.time);
M=zeros(nt,3);
se=M;
summary=zeros(3,8);
for k=1:3
    r=runs{k};
    M(:,k)=mean(r.energy,2);
    se(:,k)=std(r.energy,0,2)/sqrt(p.samples);
    ds=r.integral-proposed.integral;
    summary(k,:)=[mean(r.terminal),std(r.terminal)/sqrt(p.samples),mean(r.integral), ...
        std(r.integral)/sqrt(p.samples),mean(r.impulses),mean(r.cancellations), ...
        mean(ds),std(ds)/sqrt(p.samples)];
end
% A shared-path, doubled-resolution check for the three ablations.
% Each pair of conditional half increments sums to the coarse increment.
pCheck=p;
pCheck.samples=64;
pFine=pCheck;
pFine.dt=p.dt/2;
finePaths.switches=paths.switches(1:pCheck.samples);
coarseD=paths.dB(:,1:pCheck.samples);
rsBridge=RandStream('mt19937ar','Seed',p.seed+2);
firstHalf=0.5*coarseD+sqrt(p.dt/4)*randn(rsBridge,size(coarseD));
finePaths.dB=zeros(2*size(coarseD,1),size(coarseD,2));
finePaths.dB(1:2:end,:)=firstHalf;
finePaths.dB(2:2:end,:)=coarseD-firstHalf;
refinement=zeros(3,4);
for k=1:3
    f=simulate_bounded_wait(pFine,finePaths,names{k},NaN);
    r=runs{k};
    Iold=mean(r.integral(1:pCheck.samples)); Inew=mean(f.integral);
    refinement(k,:)=[Iold,Inew,abs(Inew/Iold-1),mean(f.impulses-r.impulses(1:pCheck.samples))];
end
assert(all(refinement(:,3)<0.02));
T=array2table(summary,'VariableNames',{'MeanTerminal','SETerminal','MeanIntegrated','SEIntegrated', ...
    'MeanImpulses','MeanCancellations','MeanPairedIntegratedExcess','SEPairedIntegratedExcess'});
T=addvars(T,displayNames,'Before',1,'NewVariableNames','Method');
writetable(T,fullfile(folder,'bw_metrics.csv'));
disp(T);
C=table(["full_P";"full_I_lower";"full_I_upper";"delay_P";"delay_I";"impulse"],cert.residual, ...
    'VariableNames',{'Condition','SignedResidual'});
writetable(C,fullfile(folder,'bw_certificate.csv'));
disp(C);
disp(struct('initialGenerator',cert.initialGenerator,'initialClockGenerator',cert.initialClockGenerator, ...
    'minM',cert.minM, ...
    'maxM',cert.maxM,'Qfp',cert.Qfp,'Qfi',cert.Qfi,'boundConstant',cert.boundConstant));
ref=array2table(refinement,'VariableNames',{'Integrated_dt','Integrated_dtHalf','RelativeDifference','DifferenceInMeanEvents'});
ref=addvars(ref,displayNames,'Before',1,'NewVariableNames','Method');
writetable(ref,fullfile(folder,'bw_refinement.csv')); disp(ref);
curve=table(proposed.time,M(:,1),se(:,1),M(:,2),se(:,2),M(:,3),se(:,3), ...
    cert.boundConstant*exp(-p.lambda*proposed.time), ...
    'VariableNames',{'time','stic','SEstic','zero_floor','SEzeroFloor','noImpulses','SEnoImpulses','certificate'});
writetable(curve,fullfile(folder,'bw_mean_curves.csv'));
eventTable=array2table(proposed.firstEvents,'VariableNames',{'time','preEnergy','postEnergy'});
writetable(eventTable,fullfile(folder,'bw_first_events.csv'));
sample=table(proposed.time,proposed.energy(:,1),'VariableNames',{'time','squaredError'});
writetable(sample,fullfile(folder,'bw_first_path.csv'));
writetable(array2table(proposed.firstCandidates,'VariableNames',{'reference','candidate','sampledQ'}), ...
    fullfile(folder,'bw_first_candidates.csv'));
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
