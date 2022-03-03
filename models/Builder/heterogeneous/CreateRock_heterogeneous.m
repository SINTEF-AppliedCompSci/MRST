function model = CreateRock_heterogeneous(model, f)
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
    coreLength = model.experiment.geometry.length.value;
    poro = experiment.rock.poro.value * ones(G.cells.num, 1);
    perm = experiment.rock.perm.value * ones(G.cells.num, 1);
  
    if (isfield(simulation,'bCells'))
        if isfield(experiment.rock.poro, 'porosity_profile')
            porosity_profile_obs = experiment.rock.poro.porosity_profile;
            length_obs = porosity_profile_obs(:,1); % care for the units
            porosity_obs = porosity_profile_obs(:,2);
            x = [0;G.cells.centroids(2:end-1,1) - 2 * G.cells.centroids(1,1);coreLength];
            poro = interp1(length_obs / 100, porosity_obs, x); 
        end
        % figure; plot(x,poro); xlabel("Distance (m)"); ylabel("Porosity")
        perm = perm .* f .^2 .* poro ./ experiment.rock.poro.value;
        % figure; plot(x,perm); xlabel("Distance (m)"); ylabel("Permeability")
    end

    rock = makeRock(G, perm, poro); 
    rock.pv = poreVolume(G, rock); 
    if (isfield(simulation,'bCells'))
        rock.regions = struct('saturation', grid.satNum);
    end
    model.rock = rock;
end