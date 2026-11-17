function export_timing_tables()
% Export the manuscript timing tables from the archived simulation records.
root = fileparts(mfilename('fullpath'));
raw = readtable(fullfile(root,'outputs','example2_trigger_episode_log.csv'));
ell=20; varpi=.01; chi=10; H=1; dt=.001;
q=exp(-chi*H)/(expm1(ell*H)-varpi);
d=arrayfun(@(Q) fzero(@(s) max(Q,q)*(expm1(ell*s)-varpi)-exp(-chi*s),[0,H]),raw.Ptheta);
candidate=raw.theta+dt*ceil(d/dt);
action=repmat("B",height(raw),1); action(candidate<raw.zeta)="A";
T=table(raw.EpisodeID,raw.theta,candidate,raw.zeta,raw.Ptheta,action, ...
    'VariableNames',{'Number','Reference','Candidate','NextSwitch','SampledEnergy','Action'});
writetable(T,fullfile(root,'outputs','table2_schedule.csv'));
% Candidates beyond the finite simulation horizon remain finite in Table 2.
assert(max(abs(candidate-[.023;.055;.099;.157;.233;.624;.749;1.347;1.619;3.055;3.604]))<1e-10);
p=bw_parameters(); paths=bw_paths(p); firstSwitch=paths.switches{1}(1);
raw=readtable(fullfile(root,'bw_first_candidates.csv'));
raw=raw(raw.reference<firstSwitch,:);
action=repmat("A",height(raw),1);
action(raw.sampledQ<=p.q)="C";
action(raw.candidate>=firstSwitch)="B";
T=table((1:height(raw))',raw.reference,raw.candidate,repmat(firstSwitch,height(raw),1), ...
    raw.candidate-raw.reference,raw.sampledQ,action, ...
    'VariableNames',{'Number','Reference','Candidate','NextSwitch','WaitingTime','SampledEnergy','Action'});
writetable(T,fullfile(root,'outputs','table3_schedule.csv'));
assert(isequal(action,["A";"A";"C";"C";"C";"C";"B"]));
end
