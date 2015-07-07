% This script is a combination of existing scripts and routines in
% mrst-co2lab, for the purpose of analyzing a user-specified formation, and
% estimating its CO2 storage capacity.

% If using formations found in Barents Sea, this script assumes you have
% downloaded the datafiles and placed them under
% mrst-core/examples/data/CO2Atlas. Format of filename is as follows:
% "Name_thickness" and "Name_top", where Name is the formation name
% (however not followed by 'fm').

moduleCheck('co2lab', 'mex',  'libgeometry', 'opm_gridprocessing');


%% 1. Load and process a formation.
N = 1;
name = 'bjarmeland'; %'Sandnesfm'; %'Utsirafm'; %'Sandnesfm';  %'Johansenfm';
[grdecl] = getAtlasGrid(name, 'coarsening',N);
G = mprocessGRDECL(grdecl{1});  % use mprocessGRDECL if its big.  % requires mex module
G = mcomputeGeometry(G);        % function requires libgeometry module

[Gt, rock] = getFormationTopGrid(name, N);
rock.poro = unique(rock.poro); % a single value. Use full porosity map instead?
if isnan(rock.poro)
    fprintf('Porosity required: Update getAtlasGrid() to include porosity of %s formation.\n', name);
    return
end

% Note, one could also use the following function, however this function
% does not return rock data as output:
% Gt = topSurfaceGrid(G);


%% 2. Visualize formation grid and top surface.
figure;
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
plotCellData(Gt, Gt.cells.z) % depth is in units of meters

view(3);
set(gca,'DataAspect',[1 1 1/100]);
hcb = colorbar;
hcb.TickLabels = hcb.Ticks/1000;
hcb.Label.String = 'Depth below sealevel, km';
xlabel('x-coordinate'); ylabel('y-coordinate'); zlabel('z, m'); grid



%% 3. Analyze traps with spill-point analysis, and plot structural traps.
ta = trapAnalysis(Gt, false); % true for cell-based method

plotCellData(Gt, ta.traps, ta.traps~=0, 'FaceColor', 'r', 'EdgeColor','none')
bf = boundaryFaces(Gt, ta.traps~=0);
plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);

ta_volumes = volumesOfTraps(Gt, ta);

fprintf('Coarsening level %d:\n', N);
fprintf('  Num. global traps: %d\n', numel(ta_volumes));
fprintf('  Total trap volume: %e m3\n', sum(ta_volumes));
fprintf('  Avg. global trap size: %e m3\n', mean(ta_volumes));

title({[grdecl{1}.name ', coarsening level=' num2str(N)]; ...
    [num2str(numel(ta_volumes)) ' structural traps shown in red']}, 'FontSize', 14)

% We can use plotCellData and color the faces according to the trap volume
% in km^3. The cell of Gt is associated with a particular trap by ta.traps,
% and is 0 otherwise. The volume of a particular trap is given by
% ta_volumes.
figure;
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')

trapcells = ta.traps~=0;
cellsTrapVol = zeros(Gt.cells.num,1);
cellsTrapVol(trapcells) = ta_volumes(ta.traps(trapcells));
plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0)

set(gca,'DataAspect',[1 1 1/100])
hcb = colorbar;
hcb.TickLabels = hcb.Ticks;
hcb.Label.String = 'Trap Volume, km^3';
xlabel('x-coordinate'); ylabel('y-coordinate'); zlabel('z, m'); grid

% Note that it may be more useful to color the faces according to the
% trapped mass of CO2. This requires using the CO2 density data which can
% vary with depth (since pressure and temperature gradients exist in the
% formation), to determine the cooresponding CO2 mass trapped by the traps.
% Porosity might also be required, in order to compute both trapping
% volumes and CO2 mass trapping.
% See functions which produce Utsira plots of traps, catchment area, etc:
%   co2 = CO2props();
%   utsiraCapacity('rhoCfun', @co2.rho);
% The above functions also compute breakdown of CO2 trapping types
% (structural, residual, and dissolution) without VE simulation.


%% 3. b) get trapping breakdown (structural, residual, dissoluion)
% The following section computes the breakdown of structural, residual, and
% dissolution trapping in the formation, assuming CO2 is injected and fills
% the formation completely, i.e., all spill paths are utilized. The
% contents of this section closely follows the implementation found in
% exploreCapacity(), which is a GUI. Here, the ranges of the input
% parameters (such as temperature gradient, etc...) are used to compute the
% range of trapping capacity (such as total and structural). Note: some
% fluid parameters are fixed to those used in exploreCapacity.


% Range of values:
tgrad_range             = 10:1:50; % degrees C (not Kelvin)
press_deviation_range   = -50:2:50;
res_co2_sat_range       = 0:0.01:0.5;
res_brine_sat_range     = 0:0.01:0.5;
diss_max_range          = [0:2:100];

% computes breakdown given range of one parameter (other parameters will be
% set to their default values)
[out_tgrad]     = exploreParameterRanges(Gt, rock, ta, 'tgrad',           tgrad_range);
[out_press]     = exploreParameterRanges(Gt, rock, ta, 'press_deviation', press_deviation_range);
[out_resco2]    = exploreParameterRanges(Gt, rock, ta, 'res_co2_sat',     res_co2_sat_range);
[out_resbri]    = exploreParameterRanges(Gt, rock, ta, 'res_brine_sat',   res_brine_sat_range);
[out_dismax]    = exploreParameterRanges(Gt, rock, ta, 'dis_max',         diss_max_range);



% visualization of trends
% TODO: make following implementation more concise (i.e., put in loop)
figure; set(gcf,'Position',[100 600 2500 400])
subplot(1,5,1)
plot(tgrad_range, out_tgrad.tot_trap_capa_sum, 'bx', ...
    tgrad_range, out_tgrad.strap_mass_co2_sum, 'ro', ...
    tgrad_range, out_tgrad.btrap_mass_co2_res_sum, 'g+', ...
    tgrad_range, out_tgrad.btrap_mass_co2_dis_sum, 'mv');
legend('total','structural','residual','dissolution')
xlabel({'Temp Gradient, degrees per kilometer';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,2)
plot(press_deviation_range, out_press.tot_trap_capa_sum, 'bx', ...
    press_deviation_range, out_press.strap_mass_co2_sum, 'ro', ...
    press_deviation_range, out_press.btrap_mass_co2_res_sum, 'g+', ...
    press_deviation_range, out_press.btrap_mass_co2_dis_sum, 'mv');
legend('total','structural','residual','dissolution')
xlabel({'Pressure deviation from hydrostatic, percent';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,3)
plot(res_co2_sat_range, out_resco2.tot_trap_capa_sum, 'bx', ...
    res_co2_sat_range, out_resco2.strap_mass_co2_sum, 'ro', ...
    res_co2_sat_range, out_resco2.btrap_mass_co2_res_sum, 'g+', ...
    res_co2_sat_range, out_resco2.btrap_mass_co2_dis_sum, 'mv');
legend('total','structural','residual','dissolution')
xlabel({'Residual CO2 saturation';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)

subplot(1,5,4)
plot(res_brine_sat_range, out_resbri.tot_trap_capa_sum, 'bx', ...
    res_brine_sat_range, out_resbri.strap_mass_co2_sum, 'ro', ...
    res_brine_sat_range, out_resbri.btrap_mass_co2_res_sum, 'g+', ...
    res_brine_sat_range, out_resbri.btrap_mass_co2_dis_sum, 'mv');
legend('total','structural','residual','dissolution')
xlabel({'Residual Brine saturation';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,5)
plot(diss_max_range, out_dismax.tot_trap_capa_sum, 'bx', ...
    diss_max_range, out_dismax.strap_mass_co2_sum, 'ro', ...
    diss_max_range, out_dismax.btrap_mass_co2_res_sum, 'g+', ...
    diss_max_range, out_dismax.btrap_mass_co2_dis_sum, 'mv');
legend('total','structural','residual','dissolution')
xlabel({'Maximum Dissolution';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)


% visualization of trends (only structural)
% TODO: make following implementation more concise (i.e., put in loop)
figure;  set(gcf,'Position',[100 100 2500 400])
subplot(1,5,1)
plot(tgrad_range, out_tgrad.strap_mass_co2_sum, 'ro');
legend('structural')
xlabel({'Temp Gradient, degrees per kilometer';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,2)
plot(press_deviation_range, out_press.strap_mass_co2_sum, 'ro');
legend('structural')
xlabel({'Pressure deviation from hydrostatic, percent';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,3)
plot(res_co2_sat_range, out_resco2.strap_mass_co2_sum, 'ro');
legend('structural')
xlabel({'Residual CO2 saturation';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)

subplot(1,5,4)
plot(res_brine_sat_range, out_resbri.strap_mass_co2_sum, 'ro');
legend('structural')
xlabel({'Residual Brine saturation';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)
    
subplot(1,5,5)
plot(diss_max_range, out_dismax.strap_mass_co2_sum, 'ro');
legend('structural')
xlabel({'Maximum Dissolution';'(other parameters are default values)'});
ylabel('Trapping Capacity, Gtons');
title(name)


                



%% 4. Use method to place wells:
% a) manual
% b) use optimization (objective function) to determine well placement of N
% wells.
% c) use array of wells




%% 5. Given placed wells, compute their injection rates:
% a) compute rate such that an injection period of X years will fill the
% volume of traps along the spill-paths over a period of Y years.
% (user-specified X and Y)




%% *************************************************************************
%% 6. Run a VE simulation, and obtain CO2 inventory breakdown:
% a) compare with and without dissolution
% b) compare bdry types (open, closed, semi)




%% 7. Use objective function to get optimized injection rates at wells.


