%% Two-Phase Flow in Inclined Gravity Column
% In this example, we simulate the injection of a light fluid (CO2) into
% a heavier fluid (brine) inside an inclined sandbox.

mrstModule add incomp


%% Set up model
% To get an inclined reservoir, we manipulate the gravity direction. Since
% gravity is a persistent and global variable, it is important that we
% reset the gravity at the end of the script
n      = [ 20,  2,  100];
box_sz = [100, 10, 200];
G      = cartGrid(n, box_sz);
G      = computeGeometry(G);
rock   = makeRock(G, 0.1*darcy, 0.3);
T      = computeTrans(G, rock, 'verbose', true);
           
fluid  = initSimpleFluid('mu' , [  0.307,   0.049] .* centi*poise     , ...
                         'rho', [973    , 617    ] .* kilogram/meter^3, ...
                         'n'  , [  2    ,   2    ]);

% Redefine gravity direction
rot   = @(theta) makehgtform('xrotate',  theta(1), ...
                             'yrotate', -theta(2), ...
                             'zrotate', -theta(3));
mul   = @(a,b,n) a(1:n,1:n) * reshape(b(1:n), [], 1);
angle = [0, pi/4, 0];
gravity reset on
gravity(mul(rot(angle), gravity(), 3));

%% Solve flow problem
% Put region of CO2 at bottom of reservoir.
xr = initResSol(G, 1*barsa, 1);
d  = gravity() ./ norm(gravity);
c  = G.cells.centroids * d.' > 110;
xr.s(c) = 0;
xr = incompTPFA(xr, G, T, fluid);

cla reset
h = plotCellData(G, convertTo(xr.pressure(1:G.cells.num), barsa), ...
                 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.625);
rotate(h,[0 1 0],pi/4*90);
view([0,0]), axis tight off
colorbar, colormap(jet)

%% Prepare to plot during simulation
clf
h = plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
rotate(h,[0 1 0],pi/4*90);
view([0,0]), axis tight off
hs = plotCellData(G, xr.s, xr.s < .995, 'EdgeColor', 'none');
rotate(hs,[0 1 0],pi/4*90);
s  = linspace(0, 1, 128).';
cm = [1-s.^(13/16), 1-s.^6, s.^6];
caxis([0 1]); colormap(cm), colorbar

%% Run simulation
% For accuracy, the time step is gradually ramped up
dT = [1, 2, 2, 5, 5, 10, 15, 20, 40, 50, 50, ...
      100, 100, 200, 200, 300, 400, 500] .* day();
dT = [dT, [2, 2, 2, 4, 5, 5, 10, 10, repmat(15, [1, 54])].*year()]/200;
t = 0;
for k = 1 : numel(dT),
   xr = implicitTransport(xr, G, dT(k), rock, fluid, 'Verbose', false);

   % Check for inconsistent saturations
   assert (max(xr.s) < 1+eps && min(xr.s) > -eps);

   % Increase time and plot saturation
   t = t + dT(k);
   delete(hs)
   hs = plotCellData(G, xr.s, xr.s <.995, 'EdgeColor', 'none');
   rotate(hs,[0 1 0],pi/4*90);
   drawnow

   % Compute new flow field.
   xr = incompTPFA(xr, G, T, fluid);
end

%% NB: RESET GRAVITY
% Gravity is defined as a persistent and global variable and we therefore
% need to reset it to avoid messing with other examples
gravity reset on