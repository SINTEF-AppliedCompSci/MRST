%% Demonstrate Interactive Plotting of Fluid Properties for AD-solvers
% The AD-OO framework can interactively visualize the fluid model of a
% ReservoirModel instance. Once active, the user can interactively explore
% the different fluid properties (viscosities, relative permeabilities,
% densities) as functions of saturation and pressure. 
%
% Load the blackoil module and others to create test data.
mrstModule add ad-blackoil ad-core ad-props deckformat example-suite

%% Inspect the SPE1 fluid model
% The SPE1 fluid model is a three-phase blackoil model with solution gas.
% For more information, as well as a simulation example, see the
% <matlab:edit('blackoilTutorialSPE1.m') SPE1 blackoil tutorial>.
%
% We use a setup routine to get a ReservoirModel subclass, and pass it to
% inspectFluidModel.
%
% Note that in this specific case, the tables for undersaturated values can
% lead to negative or discontinuous values for certain high pressure
% values. This is not always easy to see directly from the tables, but by
% using the fluid inspector it is straightforward.
[G, rock, fluid, deck] = setupSPE1();
spe1 = selectModelFromDeck(G, rock, fluid, deck);

inspectFluidModel(spe1, 'pressureRange', (0:10:500)*barsa)

%% Inspect the SPE9 fluid model
% Another standard black-oil test case is the SPE9 model. For more
% information about this test case, see the 
% <matlab:edit('blackoilTutorialSPE9.m') SPE9 black-oil tutorial>.
%
% Of particular interest in this case is the non-zero capillary pressure,
% and the highly irregular relative permeability curves.
[G, rock, fluid, deck] = setupSPE9();
spe9 = selectModelFromDeck(G, rock, fluid, deck);

inspectFluidModel(spe9)

%% Set up a two-phase oil-water fluid and inspect it
% We can use the inspection utility to get a better understanding of how
% the fluid model changes when we adjust parameters. In this case, we set
% up a simple two-phase, oil-water model and add in other relative
% permeability curves with different residual values and Corey exponents.
%
% Feel free to modify the parameters and look at how the values change.
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'rho', [1000, 700], ...
                           'cR',  1e-8/barsa, ...
                           'c',   [0, 1e-4/barsa]);
srw = 0.2;
sro = 0.3;

% Fluid relative permeabilities (use name convention from SWOF keyword)
fluid.krW  = coreyPhaseRelpermAD(2, srw, 1, srw + sro);
fluid.krOW = coreyPhaseRelpermAD(4, sro, 1, srw + sro);

model_ow = TwoPhaseOilWaterModel([], [], fluid);

% Inspect model, and specify pressure range of interest
inspectFluidModel(model_ow, 'pressureRange', (250:10:500)*barsa);

%% Copyright notice

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
