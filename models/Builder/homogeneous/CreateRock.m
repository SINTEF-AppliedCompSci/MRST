function model = CreateRock(model)
% <keywords>
%
% Purpose : assign porosity and absolute pearmbility to the grid 
%
% Syntax :
%   model = CreateRock(model)
%
% Input Parameters :
%   model: struct output from creategrid function
%
% Return Parameters :
%   model: struct containing porosity and absolute permeability
%   
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
    experiment = model.experiment;
    simulation = model.simulation;
    grid       = model.grid;
    G          = grid.G;
    poro = experiment.rock.poro.value * ones(G.cells.num, 1);
    perm = experiment.rock.perm.value * ones(G.cells.num, 1);
    
    rock = makeRock(G, perm, poro); 
    rock.pv = poreVolume(G, rock); 
    if (isfield(simulation,'bCells'))
        rock.regions = struct('saturation', grid.satNum);
    end
    model.rock = rock;
end