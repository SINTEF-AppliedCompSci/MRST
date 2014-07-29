clear

g = computeGeometry(cartGrid([5, 1]));
rock = struct('perm', ones([g.cells.num, 1]), ...
              'poro', ones([g.cells.num, 1]));

gravity([1, 0]);
gravity on

fluid  = initSimpleCompFluid('mu', [1, 1], 'n', [1, 1], 'c', [1e-3, 1e-3], ...
                            'rho_ref', [2, 1], 'press_ref', 10);
Trans  = 1 ./ accumarray(g.cells.faces(:,1), ...
                         1 ./ computeTrans(g, rock), [g.faces.num, 1]);

forces = {'bc' [], 'src', [], 'wells', []};

x      = initResSolComp(g, [], fluid, 10, [2, 1]);
porvol = poreVolume(g, rock);
T      = 0;
TSTEP  = ones(1000,1)/10;
report = [];
for k = 1 : numel(TSTEP),
   %load('previous_state')
   DT = TSTEP(k);
   [x,dt, report] = impesTPFA(x, g, Trans, fluid, DT, porvol, forces{:}, ...
                 'ATol', 5.0e-5, 'RTol', 5.0e-13, ...
                 'EstimateTimeStep', false,'DynamicMobility',false, 'report', report);


   figure(1);
   clf;
   subplot(2,1,1);
   plotCellData(g, x.pressure);colorbar;axis equal tight;

   subplot(2,1,2);
   plotCellData(g, x.z(:,1));colorbar;axis equal tight;


   %% Update for next time step
   T = T + DT;
end
