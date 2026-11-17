function reproduce_all()
% Recompute all simulation data and regenerate Figures 2--7.
root = fileparts(mfilename('fullpath'));
addpath(root);
check_example1_unified();
redraw_early_figures(fullfile(root,'figures'),root);
run_bounded_wait();
export_timing_tables();
run_comparison();
check_growth_parameters();
end
