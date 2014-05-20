%% Simulation parameters/variables for simulation of Scenario 2 in the
%% compressibility paper.

% This is the scenario of a sloping aquifer

%% Following variables must be defined in this script:
% Gt
% rock
% tinfo  = {ref_temp,  ref_depth, temp_grad}
% mu =  [muCO2, muBrine]
% slope, slopedir
% schedule

[Gt, rock] = makeTopSurfaceGrid([200, 1, 1],      ...  % # cells
                                [80000, 3000, 150], ...  % phys. dim.
                                1000, 0.1,         ...  % depth, porosity
                                1400 * milli*darcy);    % permeability


ref_temp  = 273.15 + 8; % degrees kelvin
ref_depth = 0;          % surface used as temperature reference depth
temp_grad = 45;         % degrees per kilometer
tinfo     = {ref_temp, ref_depth, temp_grad}; 
rhoW      = 1100 * kilogram / meter^3; % density of brine
mu        = [5.36108e-5, 6.5e-4];
slope     = (1/2)*pi/180;
slopedir  = [1 0];

tnum     = 50; %60; % total number of timesteps
inum     = 0;%20; % number of injection steps
tot_time = 50 * year;

% Negative values of height will be interpreted as CO2 mass, and height will
% be adjusted accordingly
h0 = zeros(200, 1);
h0(20 : 40) = -linspace(0, 100, 21) * 4e7;
h0(40 : 60) = -linspace(100, 0, 21) * 4e7;

schedule = struct('W', [], ...
                  'step', struct('val'    , diff(linspace(0, tot_time, tnum+1)), ...
                                 'control', zeros(tnum,1)));
