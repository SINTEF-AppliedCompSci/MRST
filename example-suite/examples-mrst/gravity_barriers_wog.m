function  [description, state0, model, schedule, options, plotOptions] = gravity_barriers_wog(varargin)

description = ['Gravity segregation example with water, oil and gas, ' , ...
               'and sealing, horizontal barriers, inspired from '      , ...
               'Hamon et al, 2016, doi: 10.1016/j.cma.2016.08.009'     ];
if nargout == 1, return; end

options = struct('nstep'    , 100  , ...
             'useRampup', true, ...
             'partitionType', 'voronoi');
options = merge_options(options, varargin{:});

gravity reset on
dims  = [100, 1, 100];
pdims = [100, 1, 100]*meter;
G     = cartGrid(dims, pdims);
G     = computeGeometry(G);
rock  = makeRock(G, 1*darcy, 0.3);

fluid = initSimpleADIFluid('mu' , [1, 5, 1]*centi*poise, ...
                           'rho', [1500, 1000, 500], ...
                           'cr' , 1e-6/barsa, ...
                           'n'  , [2, 2, 2]);
model = GenericBlackOilModel(G, rock, fluid);

% Estimated from "Implicit Hybrid Upwind scheme for coupled multiphase
% flow and transport with buoyancy" by Hamon et al, 2016
impermeable = false(G.faces.num, 1);
[ii, jj, kk] = gridLogicalIndices(G);

xmin = min(G.nodes.coords(:, 1));
xmax = max(G.nodes.coords(:, 1));


kmax = max(kk);
kmin = min(kk);

makex = @(frac) frac*(xmax - xmin) + xmin;
makek = @(frac) kmin + ceil(frac*(kmax - kmin));

dpth = 0.16;
a =[0.05, 0.2, dpth; ...
    0.23, 0.26, dpth; ...
    0.4, 0.6, dpth; ...
    0.7, 0.8, dpth;...
    0.82, inf, dpth];

dpth = 0.23;
b = [0.18, 0.45, dpth; ...
    0.62, 0.85, dpth];

dpth = 0.38;
c = [0, 0.1, dpth; ...
     0.18, 0.26, dpth; ...
     0.38, 0.5, dpth; ...
     0.55, 0.75, dpth; ...
     0.9, 0.95, dpth];
dpth = 0.45;
d = [0.23, 0.6, dpth; ...
    0.63, 0.71, dpth];

dpth = 0.58;
e = [0.02, 0.18, dpth; ...
    0.22, 0.3, dpth; ...
    0.48, 0.6, dpth; ...
    0.7, inf, dpth];

dpth = 0.74;
f = [-inf, 0.16, dpth; ...
     0.19, 0.23, dpth; ...
     0.24, 0.54, dpth; ...
     0.58, 0.71, dpth;...
     0.75, 0.77, dpth; ...
     0.84, 0.93, dpth];

 dpth = 0.82;
 g = [0.05, 0.19, dpth; ...
      0.24, 0.43, dpth; ...
      0.58, 0.95, dpth];

ranges = [a; b; c; d; e; f; g];
rng(0);
for i = 1:size(ranges, 1)
    k = makek(ranges(i, 3));
    faces = addSealingFaces(G, 'x_range', makex(ranges(i, 1:2)), 'k_range', [k-1, k]);
    impermeable(faces) = true;
end

model.operators.T_all(impermeable) = 0;
model.operators.T = model.operators.T_all(model.operators.internalConn);

state0 = initResSol(G, 1*atm, [0, 1, 0]);

z = G.cells.centroids(:, 3);

factor = 0.1;

low  = z > (1-factor)*pdims(end);
high = z < factor*pdims(end);

state0.s(low,  :) = repmat([0, 0, 1], sum(low), 1);
state0.s(high, :) = repmat([1, 0, 0], sum(high), 1);

time = 5*year;
dt = rampupTimesteps(time, time/options.nstep, 5*options.useRampup);

schedule = simpleSchedule(dt);

% Plotting
plotOptions = {'PlotBoxAspectRatio', [1,1,1]       , ...
               'Projection'        , 'orthographic', ...
               'View'              , [0, 0]        };

end