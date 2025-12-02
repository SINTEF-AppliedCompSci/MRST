function model = SetupCentrifuge(model)
%
% DESCRIPTION: sets the gravitational forces required for the centrifuge
% experiments
%
% SYNOPSIS:
%   model = SetupCentrifuge(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - rock: rock properties like porosity and absolute permeability
%   - grid
%   - fluid: fluid properties like viscosities and densities
%   - g: gravity vector
%   - simulation: time stepping and grid cells information
%   - twoPhaseOilWaterModel: from MRST base functions
%   - experiment: core geometry and dimensions
%
% RETURNS:
%   modifies the gravity parameters in the twoPhaseOilWaterModel struct for
%   the centrifuge experiments
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
G          = model.grid.G;
rock       = model.rock;
fluid      = model.fluid;
g          = model.gravity;
simulation = model.simulation;
nCells     = G.cells.num;

if (isfield(simulation,'bCells'))
    coreLength = model.experiment.geometry.length.value;
else
    coreLength = model.experiment.geometry.length.value * (1 - 2 / nCells);
end

centRad = model.experiment.schedule.centRad.value; 
highgrav = g * ((centRad + coreLength/2)/centRad);
lowgrav  = g * ((centRad - coreLength/2)/centRad);
dgz  = linspace(lowgrav,highgrav,nCells - 1)';
dgz(1) = 0; dgz(end) = 0; 

gravity('on','x', g);
model.twoPhaseOilWaterModel = TwoPhaseOilWaterModel(G, rock, fluid);
model.twoPhaseOilWaterModel.operators.gdz = zeros(nCells-1,1);
Grad = model.twoPhaseOilWaterModel.operators.Grad(model.twoPhaseOilWaterModel.G.cells.centroids);
model.twoPhaseOilWaterModel.operators.gdz(1:end) = Grad(:,1) .* dgz;