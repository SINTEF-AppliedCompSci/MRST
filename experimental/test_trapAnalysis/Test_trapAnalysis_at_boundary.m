%% Test and correct cell-based trap analysis at boundaries

% construct this grid:
Gt = dipped_perturbed_grid('Lx', 10000, 'Ly', 5000, 'H', 50);

% load this co2 saturation profile:
% (this state was simulated without any co2 or water residuals)
load('final_state.mat')


% compute co2 heights:
[h, h_max] = computePlumeHeight(Gt, final_state, 0, 0);


% notice when node-based trap analysis is performed, the co2 heights are
% larger than the trap heights.
ta_node = trapAnalysis(Gt, false);



% notice when cell-based trap analysis is performed, the co2 heights are
% equal to the trap heights, EXCEPT at the boundary.
ta_cell = trapAnalysis(Gt, true);