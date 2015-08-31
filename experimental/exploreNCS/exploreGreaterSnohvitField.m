%% explore Greater Snohvit Field
% located in the Hammerfest Basin Aquifer, comprised of Tubaen, Nordmela
% (low permeability), and Sto formations.

mrstModule add libgeometry opm_gridprocessing


N = 5;

%% Get datasets and construct Grids (G, Gt), rock properties
% First, get grdecl with PORO, PERMX, PERMY, PERMZ and NTG data
grdecl_st = addPermPoroNtgData2grdecl( 'Stofm', N );
grdecl_nd = addPermPoroNtgData2grdecl( 'Nordmelafm', N );
grdecl_tu = addPermPoroNtgData2grdecl( 'Tubaenfm', N );
% note: perm is in mD


% Next, we process the grids (see also: processgrid()) and compute geometry
G_st = processGRDECL(grdecl_st);
G_st = mcomputeGeometry(G_st);

G_nd = processGRDECL(grdecl_nd);
G_nd = mcomputeGeometry(G_nd);

G_tu = processGRDECL(grdecl_tu);
G_tu = mcomputeGeometry(G_tu);


% We convert PORO, PERMX, etc., and NTG into a column array, and perform
% unit conversion such that perm is in m2
rock_st = grdecl2Rock(grdecl_st, G_st.cells.indexMap);
rock_st.perm = convertFrom(rock_st.perm, milli*darcy);

rock_nd = grdecl2Rock(grdecl_nd, G_nd.cells.indexMap);
rock_nd.perm = convertFrom(rock_nd.perm, milli*darcy);

rock_tu = grdecl2Rock(grdecl_tu, G_tu.cells.indexMap);
rock_tu.perm = convertFrom(rock_tu.perm, milli*darcy);
% caution: some Inf values might be present


% Finally, we construct top-surface grids (check: do tag need to be added,
% which are needed by topSurfaceGrid? -->
% G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];  ???)
[Gt_st, G_st] = topSurfaceGrid(G_st);
[Gt_nd, G_nd] = topSurfaceGrid(G_nd);
[Gt_tu, G_tu] = topSurfaceGrid(G_tu);


% Even though the rock properties for Sto, Nordmela, and Tubaen are already
% 2D (no discretization in vertical direction), we still use averageRock to
% get rock2D which contains one column array of perm, (i.e., the perm in
% the x-direction), as VE models do not handle anisotropy.
rock2D_st  = averageRock(rock_st, Gt_st);
rock2D_nd  = averageRock(rock_nd, Gt_nd);
rock2D_tu  = averageRock(rock_tu, Gt_tu);



%% Visualization to assess grids and rock properties.
figure; set(gcf,'Position',[1 1 1000 800])
% top surface grid of Sto
%plotGrid(Gt_st,'FaceColor','none','EdgeColor','y')
% three formation grids, colored differently
plotGrid(G_st,'FaceColor','r')%,'EdgeColor','none')
plotGrid(G_nd,'FaceColor','b')%,'EdgeColor','none')
plotGrid(G_tu,'FaceColor','g')%,'EdgeColor','none')
%plotCellData(Gt_st, Gt_st.cells.z,'EdgeColor','none')
view([-100 55]); grid; legend('Sto','Nordmela','Tubaen','Location','NorthEast')
%axis equal

% Get color bar limits of the data, for consistency in plotting of three
% formations. CAUTION: Inf (or NaN) should be ignored.
ntg_max = max([ max(rock2D_st.ntg(isfinite(rock2D_st.ntg))), ...
                max(rock2D_nd.ntg(isfinite(rock2D_nd.ntg))), ...
                max(rock2D_tu.ntg(isfinite(rock2D_tu.ntg)))  ]);
            
poro_max = max([ max(rock2D_st.poro(isfinite(rock2D_st.poro))), ...
                 max(rock2D_nd.poro(isfinite(rock2D_nd.poro))), ...
                 max(rock2D_tu.poro(isfinite(rock2D_tu.poro)))  ]);
            
perm_max = max([ max(rock2D_st.perm(isfinite(rock2D_st.perm))), ...
                 max(rock2D_nd.perm(isfinite(rock2D_nd.perm))), ...
                 max(rock2D_tu.perm(isfinite(rock2D_tu.perm)))  ]);
             
ntg_min = min([ min(rock2D_st.ntg(isfinite(rock2D_st.ntg))), ...
                min(rock2D_nd.ntg(isfinite(rock2D_nd.ntg))), ...
                min(rock2D_tu.ntg(isfinite(rock2D_tu.ntg)))  ]);
            
poro_min = min([ min(rock2D_st.poro(isfinite(rock2D_st.poro))), ...
                 min(rock2D_nd.poro(isfinite(rock2D_nd.poro))), ...
                 min(rock2D_tu.poro(isfinite(rock2D_tu.poro)))  ]);
            
perm_min = min([ min(rock2D_st.perm(isfinite(rock2D_st.perm))), ...
                 min(rock2D_nd.perm(isfinite(rock2D_nd.perm))), ...
                 min(rock2D_tu.perm(isfinite(rock2D_tu.perm)))  ]);
             
depth_max = max([ max(Gt_st.cells.z), max(Gt_nd.cells.z), max(Gt_tu.cells.z) ]);
depth_min = min([ min(Gt_st.cells.z), min(Gt_nd.cells.z), min(Gt_tu.cells.z) ]);

% The following figure can be compared with the plots shown in Compiled CO2
% Atlas, chp 6, pg 129.
figure;

% Sto formation
subplot(3,4,1)
plotCellData(Gt_st, Gt_st.cells.z, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([depth_min depth_max])
title('Top surface depth, meter'); ylabel('Sto formation')
subplot(3,4,2)
plotCellData(Gt_st, rock2D_st.ntg, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([ntg_min ntg_max])
title('net-to-gross')
subplot(3,4,3)
plotCellData(Gt_st, rock2D_st.poro, 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([poro_min poro_max])
title('porosity')
subplot(3,4,4)
plotCellData(Gt_st, rock2D_st.perm./(milli*darcy), 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([perm_min perm_max]./(milli*darcy))
title('permeability, mD')

% Nordmela
subplot(3,4,5)
plotCellData(Gt_nd, Gt_nd.cells.z, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([depth_min depth_max])
ylabel('Nordmela formation')
subplot(3,4,6)
plotCellData(Gt_nd, rock2D_nd.ntg, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([ntg_min ntg_max])
subplot(3,4,7)
plotCellData(Gt_nd, rock2D_nd.poro, 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([poro_min poro_max])
subplot(3,4,8)
plotCellData(Gt_nd, rock2D_nd.perm./(milli*darcy), 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([perm_min perm_max]./(milli*darcy))

% Tubaen
subplot(3,4,9)
plotCellData(Gt_tu, Gt_tu.cells.z, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([depth_min depth_max])
ylabel('Tubaen formation')
subplot(3,4,10)
plotCellData(Gt_tu, rock2D_tu.ntg, 'EdgeColor','none');  view(2); axis equal tight; colorbar; caxis([ntg_min ntg_max])
subplot(3,4,11)
plotCellData(Gt_tu, rock2D_tu.poro, 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([poro_min poro_max])
subplot(3,4,12)
plotCellData(Gt_tu, rock2D_tu.perm./(milli*darcy), 'EdgeColor','none'); view(2); axis equal tight; colorbar; caxis([perm_min perm_max]./(milli*darcy))


%% Do some trapping analysis
% We will consider the top layer of the Hammerfest Basin aquifer first,
% which is the Sto formation. Q: which trapping estimates do porosity and
% permeability impact? Where has an averaged value (for entire formation)
% been used in estimates? Now that porosity and permeability maps are
% available, how will they impact trapping estimations?

% Determine cell index of Snohvit injection location
% Physical coordinate (approx, see chp 6 pg 135 in Atlas):
wellXcoord      = 9.202e5;
wellYcoord      = 7.988e6;
dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [wellXcoord, wellYcoord]);
[v,i] = min(sum(dv.^2, 2));
wellCellIndex = i; % or Gt.cells.indexMap(i);
%[wellInd_i, wellInd_j] = ind2sub(Gt_st.cartDims, wellCellIndex);

interactiveTrapping(Gt_st, 'injpt', wellCellIndex)
view([-100 55]);

% Note that injection wells could be placed in other locations of the Sto
% formation, to exploit the full capacity of the structural traps. Also, it
% is quite clear by the outline of some structural traps where some of the
% NPD's prospects are location. Trapping analysis should find these
% prospect areas and display detail such as rock volume, pore volume,
% storage capacity, etc.

% could also use:
exploreCapacity( 'default_formation',  'Stofm',     ...
                 'grid_coarsening',     N,          ...
                 'seafloor_depth',      330*meter,  ...
                 'seafloor_temp',       4,          ...
                 'temp_gradient',       40           );
             
exploreSimulation(  'default_formation',  'Stofm',     ...
                    'grid_coarsening',     N,          ...
                    'seafloor_depth',      330*meter,  ...
                    'seafloor_temp',       4,          ...
                    'temp_gradient',       40,         ...
                    'inj_time',            30 * year,  ...
                    'inj_steps',           10,         ...
                    'mig_time',            0 * year,   ...
                    'mig_steps',           0           );
                
% note: averaged perm and poro is used in exploreSimulation() since
% getFormationTopGrid() is called which returns rock2D data with averaged
% properties. 


%% Focus on Greater Snohvit Field for injection scenario
% then cut out region belonging to Greater Snohvit Field (where coordinates
% are taken from NDP figures in Chapter 6 of Compiled Atlas).
[wellSols, states, sim_report, opt, var ] = runSnohvitInjectionScenario( Gt_st, rock2D_st );

bf = boundaryFaces(Gt_st);
plotRealVsDiscreteInjLoc(Gt_st, bf, opt.wellXcoord, opt.wellYcoord, var.wellCoord_x, var.wellCoord_y);
makeSideProfilePlots_CO2heights(G_st, Gt_st, var.wellCellIndex, states, var.fluid );

% can do other analysis of results. i.e.,
plotToolbar(G_st,states)


% preliminary trapping anaylsis on G. Snohvit

%% Cut out Greater Snohvit Field from Sto formation and refine model



