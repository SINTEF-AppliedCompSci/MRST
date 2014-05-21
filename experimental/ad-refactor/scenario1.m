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

if ~square_domain
    % the domain is a single row of cells
    [Gt, rock] = makeTopSurfaceGrid([100, 1, 1],      ...  % # cells
                                    [40000, 3000, 150], ...  % phys. dim.
                                    750, 0.1,         ...  % depth, porosity
                                    400 * milli*darcy);    % permeability
else
    % the domain is square
    [Gt, rock] = makeTopSurfaceGrid([100 100, 1], ...
                                    [40000 40000 150], ...
                                    750, 0.1,...
                                    400*milli*darcy);
end

    
ref_temp  = 273.15 + 6; % degrees kelvin
ref_depth = 0;          % surface used as temperature reference depth
temp_grad = 40;         % degrees per kilometer
tinfo     = {ref_temp, ref_depth, temp_grad}; 
rhoW      = 1050 * kilogram / meter^3; % density of brine
mu        = [5.36108e-5, 5.4e-5];
h0        = zeros(Gt.cells.num, 1);
slope     = 0;
slopedir  = [1 0];

tnum     = 60; % total number of timesteps
inum     = 19; % number of injection steps
tot_time = 60 * year;
schedule = struct('W', addWell([], Gt, rock, ceil(Gt.cells.num/2), ...
                               'type'   , 'rate'                      , ...
                               'radius' , 0.3                         , ...
                               'comp_i' , [0 0 1]                     , ...
                               'val'    , 1e7 * kilo * kilogram /year , ...
                               'name'   , 'I'), ...
                  'step', struct('val'    , diff(linspace(0, tot_time, tnum+1)), ...
                                 'control', [ones(inum,1); zeros(tnum-inum, 1)]));
