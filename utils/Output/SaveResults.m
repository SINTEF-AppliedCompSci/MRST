function SaveResults(model)
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
end