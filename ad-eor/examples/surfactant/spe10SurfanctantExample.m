%% Oil-Water-Surfactant System for a Layer of the SPE10 Model
%
mrstModule add ad-core ad-blackoil ad-eor ad-props ...
    deckformat mrst-gui spe10 ad-fi

%% Use setupSPE10_AD to Fetch an Oil-Water Model
% We pick up only one layer
%
layers = 35;
[~, model, ~] = setupSPE10_AD('layers', layers);

G = model.G;
rock = model.rock;
fluid = model.fluid;

%% Modify the Fluid Properties
%
% Setup the fluid properties for our case
%
% Setup the relative permeabilities using Corey model for the fluid without
% surfactant and saturated with surfactant
%

% Relative permeabilities - without surfactant
n      = 3;   % Corey coefficient
sWcon  = 0.2; % Residual water saturation
sOres  = 0.2; % Residual oil saturation
krWres = 0.6; % Endpoint relperm for water
krOres = 0.5; % Endpoint relperm for oil

krW = coreyPhaseRelpermAD(n, sWcon, krWres, sWcon + sOres);
krOW = coreyPhaseRelpermAD(n, sOres, krOres, sWcon + sOres);

fluid.krW = krW;
fluid.krOW = krOW;
fluid.sWcon = sWcon;
fluid.sOres = sOres;

% Relative permeabilities - with surfactant
n         = 1.5;
sWconSft  = 0.05;
sOresSft  = 0.05;
krWresSft = 1;
krOresSft = 1;

krWSft = coreyPhaseRelpermAD(n, sWconSft, krWresSft, sWconSft + sOresSft);
krOWSft = coreyPhaseRelpermAD(n, sOresSft, krOresSft, sWconSft + sOresSft);
fluid.krWSft = krWSft;
fluid.krOWSft = krOWSft;

fluid.sWconSft   = sWconSft;
fluid.sOresSft   = sOresSft;

% Remaining fluid parameters
pRef = 234*barsa;                                        % Reference pressure
bW0        = 1./(1.012);                                 % Reference water formation volume factor
cW         = 4.28e-5/barsa;                              % Compressibility coefficient
fluid.bW   = @(p) bW0*exp((p - pRef)*cW);                % Water formation volume factor
fluid.muWr = 0.48*centi*poise;                           % Water viscosity at reference pressure
fluid.cmuW = 0/barsa;                                    % Viscosibility equal to zero
fluid.muW  = @(p) fluid.muWr*exp(fluid.cmuW*(p - pRef)); % Water viscosity (constant)
bO0        = 1./(1.065);                                 % Reference oil formation volume factor
cO         = 6.65e-5/barsa;                              % Oil compressibility coefficient
fluid.bO   = @(p) bO0*exp((p - pRef)*cO);                % Oil formation volume factor
fluid.muO  = @(p) 0*p + 5*centi*poise;                   % Oil viscosity (constant)
fluid.rhoWS = 962;                                       % Water density
fluid.rhoOS = 1080;                                      % Oil surface density

% Rock parameters (compressibilities)
cR = 3e-5/barsa;
fluid.cR = cR;
fluid.pvMultR = @(p)(1 + cR.*(p-pRef));


%% Setup the Surfactant Properties
% We use tabulated values. The surfactant parameters are the same as in the
% surfactant tutorials.
%


% Interfacial surface tension
% Let us define it as a given function (meant to be a rough interpolation of
% the data in surfac.inc)
fluid.ift = @(c) ((0.05*exp(-17*c) + 1e-6)*Newton/meter);

% Interpolation factor
miscfact = [[-10; -5.5;   -4; -3; 2], ...
            [  0;    0;  0.5;  1; 1] ...
           ];
miscfact = extendTab(miscfact); % extend to constant values.
fluid.miscfact = @(Nc, varargin) interpReg({miscfact}, Nc, {':'});

% Viscosity multiplier
muWSft = [[   0;  30; 100], ...
          [0.61; 0.8;   1] ...
         ];
muWSft(:, 2) = muWSft(:, 2)*centi*poise;
muWSft = extendTab(muWSft); % extend to constant values.
fluid.muWSft = @(c) interpReg({muWSft}, c, {':'});

% Adsorption function
surfads = [[0;    1;   30;  100], ...
           [0; 5e-4; 5e-4; 5e-4] ...
          ];
surfads = extendTab(surfads); % extend to constant values.
fluid.surfads = @(c) interpReg({surfads}, c, {':'});

% Adsorption index
fluid.adsInxSft = 1;

% Rock density (used to compute adsoprtion)
fluid.rhoRSft = 2650*kilo*gram/meter^3;

model = OilWaterSurfactantModel(G, rock, fluid);

%% Define the Wells
% The wells are set to operate with a high rate and no pressure limit.
% Hence, the pressure will rise to far above what can be used in real
% operational settings. However, the main point of the example is to force
% flow through large parts of the heterogeneous reservoir and observe the
% evolution of the displacement fronts.
injeIJ = [59  17];        % Location of injection well
prodIJ = [ 2 194];        % Location of production well
rate   = 2*meter^3/day;   % Injection rate
bhp    = 200*barsa;       % Pressure at production well
nz     = G.cartDims(3);

W = [];
% Set up injection well
W = verticalWell(W, G, rock, injeIJ(1), injeIJ(2), 1:nz, ...
                 'Type'   , 'rate', ...
                 'Val'    , rate, ...
                 'Radius' , 0.1, ...
                 'Comp_i' , [1, 0], ...
                 'name'   , 'INJE', ...
                 'Sign'   , 1);
% Set up production well
W = verticalWell(W, G, rock, prodIJ(1), prodIJ(2), 1:nz, ...
                 'Type'   , 'bhp', ...
                 'Val'    , 100*barsa, ...
                 'Radius' , 0.1, ...
                 'Comp_i' , [0, 1], ...
                 'name'   , 'PROD', ...
                 'Sign'   , -1);


%% Setup the Schedule
%
% We simulate the formation of a surfactant plug.
% Three periods:
% 1) water only
% 2) water + surfactant
% 3) water only

[W.surfact] = deal(0);
control(1).W = W;
[W([W.sign] > 0).surfact] = 50*kilogram/meter^3;
control(2).W = W;

surfinj_start_time = 1000*day;
surfinj_stop_time  = 1500*day;
end_time           = 3000*day;

dt = 10*day;
val1 = linspace(0, surfinj_start_time, round(surfinj_start_time/dt));
val2 = linspace(surfinj_start_time, surfinj_stop_time, round( ...
    (surfinj_stop_time - surfinj_start_time)/dt));
val3 = linspace(surfinj_stop_time, end_time, round( ...
    (end_time - surfinj_stop_time)/dt));

step.val     = [diff(val1'); ...
                diff(val2'); ...
                diff(val3')];
step.control = [  ones(numel(val1)-1, 1); ... 
                2*ones(numel(val2)-1, 1); ... 
                  ones(numel(val3)-1, 1)];
schedule.step    = step;
schedule.control = control;

schedule = refineSchedule(0, day*ones(10, 1), schedule);

%% Setup the initial state
%
state0 = initResSol(G, bhp, [sWcon, 1 - sWcon]);
state0.c      = zeros(G.cells.num, 1);
state0.cmax   = state0.c;
state0.ads    = computeEffAds(state0.c, 0, fluid);
state0.adsmax = state0.ads;

%% visualize the model properties
%
example_name = 'spe10';
vizSurfactantModel();

%% Run the simulation
%
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, ...
                      'plotReservoir', false);
[wellSols, states] = simulateScheduleAD(state0, model, schedule, 'afterStepFn', ...
                                        fn);

%% Inspect the results
%
figure, plotToolbar(G, states,'field','s:1'); plotWell(G,W,'height',.5);
view(-10,40); axis tight

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2016 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
