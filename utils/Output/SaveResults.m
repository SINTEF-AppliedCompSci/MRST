function SaveResults(model)
%
% DESCRIPTION: saves the simulation results into an output file
%
% SYNOPSIS:
%   SaveResults(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - output: settings related to the output reports
%
% RETURNS:
%   outputs the simulation results in a report .txt and excel file
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
output = model.output;
if(isfield(output,'quantities'))
    if(output.quantities.include)        
        quantities = output.quantities;
        SaveOutput(model,fullfile(quantities.filePath,quantities.fileName));            
    end
end
if(isfield(output,'satProfile'))
    if(output.satProfile.include)
        satProfile = output.satProfile;
        SaveSatProfile(model,fullfile(satProfile.filePath,satProfile.fileName));           
    end
end    
if(isfield(output,'saveConfig'))
    if(output.saveConfig.include)        
        saveConfig = output.saveConfig;
        SaveConfig(model,fullfile(saveConfig.filePath,saveConfig.fileName));
    end
end