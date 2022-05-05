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
        if (~isfield(quantities,'filePath') || ~isfolder('filePath'))
            mkdir('Results')
            quantities.filePath = strcat(pwd,'\Results\');
        end
        if ~isfield(quantities,'fileName')
            quantities.fileName = 'output.txt';
        end
        SaveOutput(model,fullfile(quantities.filePath,quantities.fileName));            
    end
end
if(isfield(output,'satProfile'))
    if(output.satProfile.include)
        satProfile = output.satProfile;
        if (~isfield(satProfile,'filePath') || ~isfolder('filePath'))
            mkdir('Results')
            satProfile.filePath = strcat(pwd,'\Results\');
        end
        if ~isfield(satProfile,'fileName')
            satProfile.fileName = 'satProfile.xlsx';
        end
        SaveSatProfile(model,fullfile(satProfile.filePath,satProfile.fileName));           
    end
end    
if(isfield(output,'saveConfig'))
    if(output.saveConfig.include)        
        saveConfig = output.saveConfig;
        if (~isfield(saveConfig,'filePath') || ~isfolder('filePath'))
            mkdir('Results')
            saveConfig.filePath = strcat(pwd,'\Results\');
        end
        if ~isfield(saveConfig,'fileName')
            saveConfig.fileName = 'config.txt';
        end
        SaveConfig(model,fullfile(saveConfig.filePath,saveConfig.fileName));
    end
end