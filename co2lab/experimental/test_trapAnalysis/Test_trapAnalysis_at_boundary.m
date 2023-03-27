%% Test and correct cell-based trap analysis at boundaries

% construct this grid:
Gt = dipped_perturbed_grid('Lx', 10000, 'Ly', 5000, 'H', 50);

% load this co2 saturation profile:
% (this state was simulated without any co2 or water residuals)
states=load('final_state.mat')
final_state=states.final_state;


% compute co2 heights:
[h, h_max] = computeHeightOfPlume(Gt, final_state, 0, 0);


% notice when node-based trap analysis is performed, the co2 heights are
% larger than the trap heights.
ta_node = trapAnalysis(Gt, false);



% notice when cell-based trap analysis is performed, the co2 heights are
% equal to the trap heights, EXCEPT at the boundary.
ta_cell = trapAnalysis(Gt, true);

%%
th_c=zeros(Gt.cells.num,1);
%th_n=zeros(Gt.cells.nu,1);
th_c(ta_cell.traps>0)=ta_cell.trap_z(ta_cell.traps(ta_cell.traps>0));
th_c(ta_cell.traps>0)=th_c(ta_cell.traps>0)-Gt.cells.z(ta_cell.traps>0);
%
figure()
plotCellData(Gt,th_c),colorbar
%%
figure()
plotCellData(Gt,h),colorbar
