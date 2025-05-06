%% Two-Phase Flow in Inclined Gravity Column
% In this example, we simulate the injection of a light fluid (CO2) into
% a heavier fluid (brine) inside an inclined sandbox.

mrstModule add incomp


%% Set up model
% To get an inclined reservoir, we manipulate the gravity direction. Since
% gravity is a persistent and global variable, it is important that we
% reset the gravity at the end of the script
theta  = 40;
height = 40;
exmpl  = 2;
n      = [ 20,  2,  100];
box_sz = [100, 10, 200];

% Grid
G  = cartGrid(n, box_sz);
G  = computeGeometry(G);
CG = cartGrid([1 1 1],box_sz); % used to create outline of sandbox

% Petrophysical data
if exmpl == 1
    rock   = makeRock(G, 0.1*darcy, 0.3);
else
    load rndseed.mat; rng(S);
    b = log(milli*darcy);
    a = (log(darcy)-b)/(.4 - .05);
    p = gaussianField(G.cartDims, [0.05 0.4], [3 1 11], 4.5);
    K = exp(a*(p-.05)+b);
    rock = makeRock(G, K(:), p(:));
end
T  = computeTrans(G, rock, 'verbose', true);

% Fluid
fluid  = initSimpleFluid('mu' , [  0.307,   0.049] .* centi*poise     , ...
                         'rho', [973    , 617    ] .* kilogram/meter^3, ...
                         'n'  , [  2    ,   2    ]);

% Redefine gravity direction
R = makehgtform('yrotate',-pi*theta/180);
gravity reset on
gravity( R(1:3,1:3)*gravity().' );

% Create special colormap
s  = linspace(0, 1, 64).';
cm = [1-s.^(13/16), 1-s.^6, s.^6];

%% Set initial data and compute pressure distribution
% Put region of CO2 at bottom of reservoir.
xr = initResSol(G, 1*barsa, 1);
d  = gravity() ./ norm(gravity);
dc = G.cells.centroids * d.';
xr.s(dc>max(dc)-height) = 0;
xr = incompTPFA(xr, G, T, fluid);

%% Plot initial data
clf
h = plotGrid(CG, 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2);
rotate(h,[0 1 0],theta);
view([0,0])
hs = plotCellData(G, xr.s, xr.s < .995, 'EdgeColor', 'none');
rotate(hs,[0 1 0],theta);
caxis([0 1]); colormap(cm), axis equal tight off

% dPlot = [20 100 250 500 1000 1500 2500 inf]*day;
% print('-dpng', '-r0', sprintf('buoy-%d-%02d.png',exmpl, 0));

%% Run simulation
% For accuracy, the time step is gradually ramped up
dT = [.5, .5, 1, 1, 1, 2, 2, 2, 5, 5, 10, 10, 15, 20, repmat(25,[1,97])].*day;
[t, ip] = deal(0,1);
for k = 1 : numel(dT)
   xr = implicitTransport(xr, G, dT(k), rock, fluid, 'Verbose', false);

   % Check for inconsistent saturations
   assert (max(xr.s) < 1+eps && min(xr.s) > -eps);

   % Increase time and plot saturation
   t = t + dT(k);
   delete(hs)
   hs = plotCellData(G, xr.s, xr.s <.995, 'EdgeColor', 'none');
   rotate(hs,[0 1 0],theta);
   title(sprintf('%.2f days', t/day));
   drawnow

   %{
   if t>=dPlot(ip)-eps,
       colorbar off; title([]);
       drawnow;
       print('-dpng', '-r0', sprintf('buoy-%d-%02d.png',exmpl, ip));
       colorbar, title(sprintf('%.2f days', t/day));
       ip = ip+1;
   end
   %}

   % Compute new flow field.
   xr = incompTPFA(xr, G, T, fluid);
end

%% NB: RESET GRAVITY
% Gravity is defined as a persistent and global variable and we therefore
% need to reset it to avoid messing with other examples
gravity reset on