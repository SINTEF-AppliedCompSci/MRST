%% SPE 6th comparative study (water injection scenario) with the dual porosity-dual permeability model
%
% Modifications with regards to the original problem formulation of
% Firoozabadi and Thomas, 1990:
%   1. Consider only oil and water
%   2. Represent oil as dead oil with a typical compressibility
%   3. Matrix permeability is increased by the factor of 10 to achieve more 
%      pronounced differences with the dual porosity-dual permeability
%      model.
%
% Reference
% A. Firoozabadi, L. K. Thomas, Sixth SPE Comparative Solution Project: 
%   Dual-Porosity Simulators, Journal of Petroleum Technology 42 (06): 710–763.
%   DOI: https://doi.org/10.2118/18741-PA.

mrstModule add ad-core ad-props ad-blackoil dual-porosity deckformat

%% Choose between the dual porosity or dual porosity-dual permeability models
model_name = 'DP';
%model_name = 'DPDP';

%% Define the geometry
[Nx, Ny, Nz] = deal(10, 1, 5);
dx = 200*ft;
dy = 1000*ft;
dz = 50*ft;
G = cartGrid([Nx Ny Nz],[Nx*dx Ny*dy Nz*dz]);
G = computeGeometry(G); 

%% Rock properties
phim = 0.29;             % Matrix porosity
% Km = 1*milli*darcy;     % Matrix permeability
Km = 10*milli*darcy;     % Increased matrix permeability
phif = 0.01;             % Fracture porosity
Kf1 = 10*milli*darcy;    % Fracture permeability for layers 1-2
Kf2 = 90*milli*darcy;    % Fracture permeability for layers 3
Kf3 = 20*milli*darcy;    % Fracture permeability for layers 4-5

Kf = ones(Nx * Ny * Nz, 1);
Kf(1:2*Nx*Ny) = Kf1;
Kf(2*Nx*Ny+1:3*Nx*Ny) = Kf2;
Kf(3*Nx*Ny+1:5*Nx*Ny) = Kf3;

% Create the rock structures for fractures and the matrix
rock_fracture = makeRock(G, Kf, phif);
rock_matrix = makeRock(G, Km, phim);

% Set rock compressibility
pref = 6000 * psia;
cr = 3.5e-6 / psia;
rock_fracture.cr = cr;
rock_fracture.pref = pref;

rock_matrix.cr = cr;
rock_matrix.pref = pref;

% plotCellData(G, rock_fracture.perm/milli/darcy), view(3), colorbar

%% Fluid properties

% Water and oil viscosities, densities, FVFs, and compressiblities
mu = deal([0.35, 0.2178] * centi*poise);
rho = deal([65, 51.14] * pound/ft^3);
fvf = deal([1.07, 1.8485]);
c = deal([3.5E-06, 1e-5] / psia);

% Water-oil relative permeabilities and capillary pressure data for the 
% fractures' and for the matrix domains, respectively. The first doman has to 
% correspond to the fractures' rock functions.
SWOF = {
[[0     0       1       0];
 [1     1       0       0]];
[[0.2   0       1       1];
 [0.25  0.005   0.86    0.5];
 [0.3   0.01    0.723	0.3];
 [0.35  0.02    0.6		0.15];
 [0.4   0.03    0.492	0.0];
 [0.45  0.045   0.392	-0.2];
 [0.5   0.06    0.304	-1.2];
 [0.6   0.11    0.154	-4.];
 [0.7   0.18    0.042	-10.];
 [0.75  0.23    0		-40.];
 [1     1		0		-100]]};

% The residual saturations for the fractures' and matrix domains
Swrf = SWOF{1}(1, 1);
ind = find(SWOF{1}(:, 3) == 0);
Snrf = 1 - SWOF{1}(ind(1), 1);

Swrm = SWOF{2}(1, 1);
ind = find(SWOF{2}(:, 3) == 0);
Snrm = 1 - SWOF{2}(ind(1), 1);

% Convert pc to Pa
for r = 1:length(SWOF)
    SWOF{r}(:, 4) = SWOF{r}(:, 4) * psia;
end

% Create the fluid model
fluid = struct();    

% Mock up the info on saturation regions for assignSWOF 
reg.sat = length(SWOF);
reg.optimize = 0;
reg.interp1d = @interpTable;
fluid = assignSWOF(fluid, SWOF, reg);

% Assign the rock compressibility to account for changing PV 
reg.pvt = 1;
reg.prange = [];
fluid = assignROCK(fluid, [pref cr], reg);

% Assign the water properties
fluid = assignPVTW(fluid, [pref fvf(1) c(1) mu(1) 0], reg);

% Assign the oil properties
fluid = assignPVCDO(fluid, [pref fvf(2) c(2) mu(2) 0], reg);

% Assign the surface densities
fluid.rhoWS = rho(1);
fluid.rhoOS = rho(2);
 
% Collapse the function cell arrays if only one region is present
fn = fieldnames(fluid);
for i = 1:numel(fn)
    f = fn{i};
    if iscell(fluid.(f)) && numel(fluid.(f)) == 1
        fluid.(f) = fluid.(f){1};
    end
end

% Assign different relative permeabilities and capillary pressures to fluds
% in the fractures' and in the matrix domains, all the other properties are
% the same
[fluid_fracture, fluid_matrix] = deal(fluid);

% Provide a dummy varargin argument as required in function calls from DualPorosityReservoirModel.m
fluid_fracture.krW = @(sw, varargin) fluid.krW{1}(sw);
fluid_fracture.krOW = @(so, varargin) fluid.krOW{1}(so);
fluid_fracture.pcOW = @(sw, varargin) fluid.pcOW{1}(sw);

fluid_matrix.krW = @(sw, varargin) fluid.krW{2}(sw);
fluid_matrix.krOW = @(so, varargin) fluid.krOW{2}(so);
fluid_matrix.pcOW = @(sw, varargin) fluid.pcOW{2}(sw);

% Set the residual saturations which can be used in the transfer function                            
fluid_fracture.swr = Swrf;
fluid_fracture.snr = Snrf;
fluid_matrix.swr = Swrm;
fluid_matrix.snr = Snrm;

%% Appy the vertical transmissibility multipliers in z-direction

% Compute the transmissibilities for the fractures and for the matrix
Tf = getFaceTransmissibility(G, rock_fracture);
Tm = getFaceTransmissibility(G, rock_matrix);

% Compute the face transmissibility multipliers in z-direction
mult = 0.1;
deck = struct('MULTZ', ones(prod(G.cartDims), 1)*mult);
m = computeTranMult(G, deck);

% Convert the multipliers from cell-based to faces-based indexing
ix   = [G.cells.faces(:,1), ...
              rldecode(G.cells.indexMap, diff(G.cells.facePos))];
[~, ia, ~] = unique(ix(:, 1));
m = m(ia);
% plotFaces(G, m < 1)

% Apply the multipliers
Tf = Tf .* m;
Tm = Tm .* m;


%% Setup the model

% The case includes gravity
gravity reset on;

% The optional inputdata parameter to the model constructor is interpreted
% as an Eclipse deck structure. We mock up the required format by introducing
% the dummy fields RUNSPEC, PROPS etc.
% The modified transmissibilities are passed to setupOperatorsTPFA by
% re-formatting the options in DualPorosityReservoirModel.
inputdata.trans = {Tf, Tm};
inputdata.RUNSPEC = true;
inputdata.PROPS = true;
inputdata.GRID = true;
inputdata.SOLUTION = true;
%inputdata = [];

if strcmp(model_name, 'DP')
    % Initialize the 2 phase dual porosity model using the same fluid model for 
    % the matrix and fracture regions   
    model = myTwoPhaseOilWaterModelDP(G, ...
                                {rock_fracture, rock_matrix}, ...  
                                {fluid_fracture, fluid_matrix}, ... 
                                'inputdata', inputdata ...
                                );
    disp('Running the dual porosity simulation..');
elseif strcmp(model_name, 'DPDP')
    % Initialize the 2 phase dual porosity-dual permeability model using the same fluid model for 
    % the matrix and fracture regions   
    model = TwoPhaseOilWaterModelDPDP(G, ...
                                {rock_fracture, rock_matrix}, ...  
                                {fluid_fracture, fluid_matrix}, ... 
                                'inputdata', inputdata ...
                                );
    disp('Running the dual porosity-dual permeability simulation..');                        
else
    disp('Model name not defined! Exiting..');
    return;
end          
                        
% The shape factors
sigma1 = 0.04/ft^2;    % Layers 1-2
sigma2 = 1/ft^2;       % Layers 3
sigma3 = 0.25/ft^2;    % Layers 4-5

sigma = ones(Nx * Ny * Nz, 1);
sigma(1:2*Nx*Ny) = sigma1;
sigma(2*Nx*Ny+1:3*Nx*Ny) = sigma2;
sigma(3*Nx*Ny+1:5*Nx*Ny) = sigma3;

% Matrix block heights for calculating the gravity drainage 
dzm1 = 25*ft;    % Layers 1-2
dzm2 = 5*ft;     % Layers 3
dzm3 = 10*ft;    % Layers 4-5

dzm = ones(Nx * Ny * Nz, 1);
dzm(1:2*Nx*Ny) = dzm1;
dzm(2*Nx*Ny+1:3*Nx*Ny) = dzm2;
dzm(3*Nx*Ny+1:5*Nx*Ny) = dzm3;

block_dimension = repmat([dx,dy,1],G.cells.num,1);
block_dimension(:, 3) = dzm;

% Initialize the transfer function 
model.transfer_model_object = myEclipseTwoPhaseTransferFunction('VariableShapeFactor', ...
   [sigma block_dimension]);

%% Initial conditions

% Initial pressure
p0 = 6000*psia;

% Set the initial pressure in the both continua to p0
clear state
state.pressure = ones(G.cells.num,1) * p0;          
state.pressure_matrix = ones(G.cells.num,1) * p0;  

% Set the initial water saturation to the corresponding Swr
state.s = repmat([Swrf 1-Swrf],G.cells.num,1);
state.sm = repmat([Swrm 1-Swrm],G.cells.num,1);


%% Set up and run the schedule

% Connect the wells to the fracture continuum
W = verticalWell([], G, rock_fracture, 1, 1, 1:Nz, ...
                 'type', 'rate', 'val', 1750*stb/day, ...
                 'Radius', 0.5*ft, ...
                 'WI', 2*centi*poise*stb/(day*psia), ...
                 'comp_i', [1, 0], ...
                 'name', 'I', 'dir', 'z');
W.lims.rate  =  1750*stb/day;             
W.lims.bhp  =  6100 * psia;

W = verticalWell(W, G, rock_fracture, 10, 1, 1:3, ...
                 'type', 'lrat', 'val', -1000*stb/day, ...
                 'Radius', 0.5*ft, ...
                 'WI', 2*centi*poise*stb/(day*psia), ...
                 'comp_i', [0.5, 0.5], ...
                 'name', 'P', 'dir', 'z'); 

W(end).lims.lrat  =  -1000*stb/day;
W(end).lims.bhp  =  1 * atm;

dt = 1*year;
n = 20;   

% Uniform time steps, spanning 20 years
schedule = simpleSchedule(ones(n, 1)*dt, 'W', W);

% Run the simulation
[wellSols, states] = simulateScheduleAD(state, model, schedule);


%% Plotting

% Get the MRST well solution
qo = zeros(length(wellSols), 1);
qw = zeros(length(wellSols), 1);
wct = zeros(length(wellSols), 1);
for i = 1:length(wellSols)
    qo(i) = - wellSols{i}(2).qOs / (stb/day);
    qw(i) = wellSols{i}(1).qWs / (stb/day);
    wct(i) = wellSols{i}(2).wcut;
end
t = cumsum(schedule.step.val) / year;

% Get the Eclipse well solution
if strcmp(model_name, 'DP')
    smry = readEclipseSummaryUnFmt('SPE6_WO_DP_MOD');
elseif strcmp(model_name, 'DPDP')
    smry = readEclipseSummaryUnFmt('SPE6_WO_DPDP_MOD');
else
    disp('Model name not defined! Exiting..');
    return;
end  
t_ref = smry.get(':+:+:+:+', 'TIME', ':');
t_ref = t_ref*day/year;
qo_ref = smry.get('FIELD', 'FOPR', ':');
wct_ref = smry.get('FIELD', 'FWCT', ':');

% Compare the well results
figure(1);
subplot(2, 1, 1)
plot(t, qo, 'r')
hold on
plot(t_ref, qo_ref, 'b')
ylabel('Oil rate [stb/day]')
title('Field oil production rate');
xlim([0 max(t)])
legend('MRST', 'ECLIPSE')

subplot(2, 1, 2)
plot(t, wct, 'r')
hold on
plot(t_ref, wct_ref, 'b')
ylabel('Water cut');
title('Water cut');
xlim([0 max(t)])
legend('MRST', 'ECLIPSE')
xlabel('Time [years]')

% Plot the evolution of 3D properties
disp(' ')
fig = figure(2);
for i = 1:length(states)    

    clf
    hp = uipanel('Parent', fig, 'BorderType', 'none', ...
                'TitlePosition', 'centertop', 'FontSize', 10); 
    set(hp, 'Title', ['At ' num2str(t(i)) ' years']);
    
    subplot(2, 2, 1, 'Parent', hp)
    plotCellData(G, states{i}.pressure / psia);
    plotWell(G, W)
    view(3), colorbar
    axis equal
    title('Fracture pressure [psia]')
    
    subplot(2, 2, 2, 'Parent', hp)
    plotCellData(G, states{i}.s(:,1));
    plotWell(G, W)
    view(3), colorbar
    axis equal
    title('Fracture water saturation')    
    
    subplot(2, 2, 3, 'Parent', hp)
    plotCellData(G, states{i}.pressure_matrix / psia);
    plotWell(G, W)
    view(3), colorbar
    axis equal
    title('Matrix pressure [psia]')
    
    subplot(2, 2, 4, 'Parent', hp)
    plotCellData(G, states{i}.sm(:,1));
    plotWell(G, W)
    view(3), colorbar
    axis equal
    title('Matrix water saturation')    
    
    drawnow

    inp = input('Press Enter to continue..');
 
end

%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

