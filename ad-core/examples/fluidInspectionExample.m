%% Demonstrate interactive plotting of fluid properties for AD-solvers
% The AD-OO framework can interactively visualize the fluid model of a
% ReservoirModel instance. Once active, the user can interactively explore
% the different fluid properties (viscosities, relative permeabilities,
% densities) as functions of saturation and pressure. 
%
% Load the blackoil module and others in order to create test data.
mrstModule add ad-blackoil ad-core ad-props deckformat
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
model_spe1 = selectModelFromDeck(G, rock, fluid, deck);

inspectFluidModel(model_spe1)
%% Inspect the SPE9 fluid model
% Another standard blackoil test case is the SPE9 model. For more
% information about this test case, see the 
% <matlab:edit('blackoilTutorialSPE9.m') SPE9 blackoil tutorial>.
%
% Of particular interest in this case is the non-zero capillary pressure,
% and the highly irregular relative permeability curves.
[G, rock, fluid, deck] = setupSPE9();
model_spe9 = selectModelFromDeck(G, rock, fluid, deck);

inspectFluidModel(model_spe9)
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
% Fluid relative permeabilities
fluid.krW = coreyPhaseRelpermAD(2, srw, 1, srw + sro);
fluid.krO = coreyPhaseRelpermAD(4, sro, 1, srw + sro);

model_ow = TwoPhaseOilWaterModel([], [], fluid);
% Inspect model, and specify pressure range of interest
inspectFluidModel(model_ow, 'pressureRange', (250:10:500)*barsa);
