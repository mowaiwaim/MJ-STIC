function redraw_all()
% Regenerate Figures 2--7 from the archived data without new simulations.
root = fileparts(mfilename('fullpath'));
addpath(root);
redraw_early_figures(fullfile(root,'figures'));
render_existing_curves();
plot_comparison_figure(readtable(fullfile(root,'outputs','comparison_trajectories.csv')), ...
    fullfile(root,'figures','figure7.pdf'));
end
