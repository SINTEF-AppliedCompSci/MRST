%% Simulation parameters/variables for simulation of Scenario 1 in the
%% compressibility paper.

%% Following variables must be defined in this script:
% Gt
% rock
% tinfo  = {ref_temp,  ref_depth, temp_grad}
% mu =  [muCO2, muBrine]
% slope (slopedir, theta)
% h0
% schedule
% refcell_ix - index to reference cell for computing constant rho (default: 1)

res = 100;

top_geom = +80 * cos([1:(res+1)]'/(res+1) * 2 * pi) + 780;
% ensuring top-geometry is exactly symmetric around center
top_geom = (top_geom + fliplr(top_geom))/2;

[Gt, rock] = makeTopSurfaceGrid([res, 1, 1],        ...  % # cells
                                [10000, 3000, 100], ...  % phys. dim.
                                top_geom, 0.2,     ...  % depth, porosity
                                1100 * milli*darcy);      % permeability

ref_temp  = 273.15 + 6; % degrees kelvin
ref_depth = 0;          % surface used as temperature reference depth
temp_grad = 40;         % degrees per kilometer
tinfo     = {ref_temp, ref_depth, temp_grad}; 
rhoW      = 1100 * kilogram / meter^3; % density of brine
mu        = [5.36108e-5, 6.5e-4];
h0        = zeros(100, 1);
slope     = 0;
slopedir  = [1 0];

refcell_ix = ceil(res/4);

tnum     = 200; %60; % total number of timesteps
inum     = 50;%20; % number of injection steps
tot_time = 200 * year;
schedule = struct('W', addWell([], Gt, rock, ceil(Gt.cells.num/2), ...
                               'type'   , 'rate'                      , ...
                               'radius' , 0.3                         , ...
                               'comp_i' , [0 0 1]                     , ...
                               'val'    , 2e6 * kilo * kilogram /year , ...
                               'name'   , 'I'), ...
                  'step', struct('val'    , diff(linspace(0, tot_time, tnum+1)), ...
                                 'control', [ones(inum,1); zeros(tnum-inum, 1)]));
