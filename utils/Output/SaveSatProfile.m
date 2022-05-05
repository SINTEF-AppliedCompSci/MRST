function SaveSatProfile(model, fullFile)
%
% DESCRIPTION: saves the saturation profile calculated with the simulation
% in an excel file
%
% SYNOPSIS:
%   SaveSatProfile(model, fullFile)
%
% PARAMETERS:
%   - model - main modeling struct with the saturation profiles saved inside
%   - fullFile - path to save the saturation profile
%
% RETURNS:
%   saturation profile saved in excel format
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
    dynamic    = model.dynamic;
    params     = dynamic.params;    
    G          = model.grid.G;
    coreLength = model.experiment.geometry.length.value;
    state0     = model.state0;
    
    if (isfield(model.simulation,'bCells'))
        model.Sw_profile = zeros(length(params.cumScheduleSteps) + 1, G.cells.num + 1);
        model.Sw_profile(2:end,1) = params.cumScheduleSteps;
        model.Sw_profile(1,2:end) = [0;G.cells.centroids(2:end-1,1)...
            - 2 * G.cells.centroids(1,1);coreLength];
        model.Sw_profile(2,2) = interp1(G.cells.centroids(2:3,1),state0.s(2:3,1),0,'linear','extrap');
        model.Sw_profile(2,end) = interp1(G.cells.centroids(end-2:end-1,1) - ...
            G.cells.centroids(1,1),state0.s(end-2:end-1,1),coreLength,'linear','extrap');
        model.Sw_profile(2,3:end-1) = state0.s(2:end-1,1)';
        for i = 1:length(dynamic.states)
            model.Sw_profile(i+2,3:end-1) = dynamic.states{i, 1}.s(2:end-1,1)';
            model.Sw_profile(i+2,2) = interp1(G.cells.centroids(2:3,1),dynamic.states{i, 1}.s(2:3,1)',0,'linear','extrap');
            model.Sw_profile(i+2,end) = interp1(G.cells.centroids(end-2:end-1,1) - ...
                            G.cells.centroids(1,1),dynamic.states{i, 1}.s(end-2:end-1,1),coreLength,'linear','extrap');
        end
    else
        model.Sw_profile = zeros(length(params.cumScheduleSteps) + 1, G.cells.num + 3);
        model.Sw_profile(2:end,1) = params.cumScheduleSteps;
        model.Sw_profile(1,2:end) = [0;G.cells.centroids(:,1);coreLength];
        model.Sw_profile(2,2) = interp1(G.cells.centroids(2:3,1),state0.s(2:3,1),0,'linear','extrap');
        model.Sw_profile(2,end) = interp1(G.cells.centroids(end-2:end-1,1) - ...
            G.cells.centroids(1,1),state0.s(end-2:end-1,1),coreLength,'linear','extrap');
        model.Sw_profile(2,3:end-1) = state0.s(:,1)';
        for i = 1:length(dynamic.states)
            model.Sw_profile(i+2,3:end-1) = dynamic.states{i, 1}.s(:,1)';
            model.Sw_profile(i+2,2) = interp1(G.cells.centroids(2:3,1),dynamic.states{i, 1}.s(2:3,1)',0,'linear','extrap');
            model.Sw_profile(i+2,end) = interp1(G.cells.centroids(end-2:end-1,1) - ...
                            G.cells.centroids(1,1),dynamic.states{i, 1}.s(end-2:end-1,1),coreLength,'linear','extrap');
        end
    end    
    writematrix(model.Sw_profile,fullFile)
    writematrix('location/time',fullFile,'Sheet',1,'Range','A1:A1') 
end