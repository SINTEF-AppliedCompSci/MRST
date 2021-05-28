%% Example Simulation of Model One From Tenth SPE CSP
mrstModule add ad-core ad-blackoil ad-props mrst-gui spe10
gravity reset on

%% Model Geometry
% Grid is a 100-by-1-by-20 Cartesian box with equally sized cells of
% dimensions 25-by-25-by-2.5 feet.
cartDims = [ 100, 1, 20 ];
physDims = cartDims .* [ 25, 25, 2.5 ]*ft;
G = computeGeometry(cartGrid(cartDims, physDims));

%% Petrophysical Properties
% Porosity is constant (=0.2) throughout the formation.  The permeability
% distribution is a correlated geostatistically generated field stored in a
% file supplied by the SPE.

rock = getSPE10_model_1_rock();

clf
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)))
view(3), axis tight, grid on

%% Define Sources and Sinks
% Model is produced from a single producer located at the far end of the
% model (I==100) constrained at a bottom-hole pressure of 95 Psi.  There is
% a single injector at the near end (I==1) providing pressure support.  The
% injector fills the reservoir with 6.97 cubic metres of gas per day.  Both
% wells have an internal diameter of 1 ft.

IJK = gridLogicalIndices(G);
I   = find(IJK{1}(:,1) == 1);
P   = find(IJK{1}(:,1) == G.cartDims(1));  clear IJK

W = addWell([], G, rock, I, 'Comp_i', [ 0, 1 ], 'Type', 'rate', ...
            'Val', 6.97*meter^3/day, 'Radius', (1/2)*ft, 'Dir', 'z', ...
            'Sign', +1, 'Name', 'I', 'refDepth', 0*ft);  clear I

W = addWell(W, G, rock, P, 'Comp_i', [ 1, 1 ], 'Type', 'bhp', ...
            'Val', 95*psia, 'Radius', (1/2)*ft, 'Dir', 'z', ...
            'Sign', -1, 'Name', 'P', 'refDepth', 0*ft);  clear P

%% Official Benchmark Relative Permeability Data
% Build a reduced ECLIPSE-style input deck that contains just enough
% information to construct relative permeability curves based on the
% official benchmark data.  In particular we use the fact that the relative
% permeability data is formatted in the same way as ECLIPSE's 'SGOF'
% keyword data.
kr_deck = getSPE10_model_1_relperm();

clf
plot(kr_deck.PROPS.SGOF{1}(:, 1),            ...
     kr_deck.PROPS.SGOF{1}(:, [2, 3]), '*-', ...
     'LineWidth', 2, 'MarkerSize', 5)
legend({'kr_g', 'kr_o'}, 'Location', 'Best')
xlabel('S_g'), title('Relative Permeability, Model I 10th CSP')

%% Fluid Properties
% The fluids in this simulation model are incompressible and immiscible
% with constant viscosities.  This means we can use MRST's special purpose
% fluid constructor |initSimpleADIFluid| to create the fluid object.  We
% will however need to use sampled relative permeability curves so we do
% not enter any relative permeability data in this call.
fluid = initSimpleADIFluid('mu'    , [  1, 0.01]*centi*poise, ...
                           'rho'   , [700, 1   ]*kilogram/meter^3, ...
                           'cR'    , 6.0e-5/barsa, ...
                           'phases', 'OG');

%%
% Replace the synthetic relative permeability curves created through
% function |initSimpleADIFluid| with the real benchmark values.
fluid_kr = assignSGOF(fluid, kr_deck.PROPS.SGOF, struct('sat', 1, ...
                                               'interp1d', @interpTable));
fluid.krG = fluid_kr.krG{1};
fluid.krO = fluid_kr.krOG{1};              clear fluid_kr

%%
% The <matlab:mrstModule('add','spe10') SPE 10 module> contains the special
% purpose function |getSPE10_model_1_fluid| that performes the above fluid
% manipulations so one would generally not do this manually.

%% Form Reservoir Model
% This is an incompressible, immiscible oil/gas system.
model = GenericBlackOilModel(G, rock, fluid, 'gravity', gravity, 'disgas', false,...
    'vapoil', false, 'water', false, 'oil', true, 'gas', true);

%% Initialise Formation
% Formation is initially filled with oil and the initial pressure at the
% top of the model is 100 Psi.
region = getInitializationRegionsBlackOil(model, 0, 'datum_pressure', 100*psia);
state0 = initStateBlackOilAD(model, region);

clf
plotCellData(G, convertTo(state0.pressure, psia))
view(3), colorbar(), axis tight, grid on
xlabel('x'), zlabel('Depth'), title('Initial Pressure Distribution [Psi]')

%%
% Ten years (3650*day), ramp-up time-steps.
timesteps = [ 0.1, 0.2, 0.3, 0.4, repmat(0.5, [1, 6]), ones([1, 6]), ...
              repmat(2, [1, 5]), repmat(5, [1, 8]), repmat(10, [1, 9]), ...
              repmat(20, [1, 5]), repmat(30, [1, 50]), ...
              repmat(50, [1, 38]) ]*day;

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);

%%
fn = getPlotAfterStep(state0, model, schedule, 'view', [50, 50], ...
                      'field', 's:1', 'wells', W);

[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, 'afterStepFn', fn);

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
