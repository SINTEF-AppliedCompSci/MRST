function state = Initialize(model)
%
% DESCRIPTION: sets the initial state of the model
%
% SYNOPSIS:
%   state = Initialize(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - simulation: time stepping and grid cells information
%   - experiment: saturation functions used for forward modeling
%   - grid
%
% RETURNS:
%   state: initial state of the model
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
Swi        = model.experiment.rock.Swi.value;
pini       = model.experiment.schedule.pini.value;
simulation = model.simulation;
G          = model.grid.G;
state      = initState(G, [], pini, [Swi 1-Swi]); % filled with phase 2 (oil)
if(isfield(simulation,'bCells'))
    firstCellSw    = simulation.bCells.firstCellSw;
    lastCellSw     = simulation.bCells.lastCellSw;
    state.s(1,1)   = firstCellSw.value;
    state.s(1,2)   = 1 - firstCellSw.value;           
    state.s(end,1) = lastCellSw.value;
    state.s(end,2) = 1 - lastCellSw.value;
end   
