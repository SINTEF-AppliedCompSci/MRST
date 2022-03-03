function model = Configure(path, fileName)
% <keywords>
%
% Purpose : read the required parameter for modeling from the settings
% files
%
% Syntax :
%   model = Configure(path, fileName)
%
% Input Parameters :
%   path: string directory of the "Data" folder
%   fileName: string name of the settings file
%
% Return Parameters :
%   model: struct containing the required parameters to start the modeling
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

    fprintf("Configuring from: %s\n", fileName);
    fprintf("Located at: %s\n", path);
    fprintf("Current directory: %s\n", pwd);

    full_path = fullfile(path,fileName);
    file  = fopen(full_path,'r');
    headers = ["GEOMETRY" "ROCK" "FLUID" "PROCESS" "SIMULATION"...
               "SCHEDULE" "OBSERVATION" "SATURATION" "PLOT" "OUTPUT"...
               "OBJECTIVE" "HISTORY" "MCMC"];
    str = fileread(full_path); % read entire file into string
    lines = strtrim(regexp( str,'(\r|\n)+','split')); % split by each line
    
    whiteSpace = 1;
    m = 1; k = 1;
    exitFlag = 0;
    foundHeader = 0;
    counter = 0;
    blockStart = [];
    while (~exitFlag)
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        counter = counter + 1;
        while (~isempty(find(strcmp(line,headers{m}),1)))
            blockStart = [blockStart,k];
            foundHeader = foundHeader + 1;
            if (foundHeader == length(headers))
                exitFlag = 1;
                break
            else
                m = m + 1;                
            end
        end
        k = k + 1;
    end
    blockStart = blockStart + 2;

    %% block GEOMETRY = ["LENGTH" "DIAMETER"]
    blockNo = 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            value = str2double(line(2));
            if(length(line) == 3), unit = line(3); end
            if(strcmp(line(1),strcat('LENGTH')))
                geometry.length.include = true;
                geometry.length.inputValue = value;
                geometry.length.inputUnit = string(unit);
                geometry.length.value = value * Convert(unit);
                break
            elseif(strcmp(line(1),strcat('DIAMETER')))
                geometry.diameter.include = true;
                geometry.diameter.inputValue = value;
                geometry.diameter.inputUnit = string(unit);
                geometry.diameter.value = value * Convert(unit);
                break
            else
                break
            end
        end
    end
    experiment.geometry = geometry;
      
    %% block ROCK = ["POROSITY" "PERMEABILITY" "INITIALWATER"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    rock.heterogeneous = false;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if length(line)>1; value = str2double(line(2)); end
            if(length(line) == 3), unit = line(3); end
            if(strcmp(line(1),strcat('POROSITY')))
                rock.poro.include = true;               
                rock.poro.inputValue = value;
                rock.poro.value = value;             
                break
            elseif(strcmp(line(1),strcat('PERMEABILITY')))
                rock.perm.include = true;               
                rock.perm.inputValue = value;
                rock.perm.inputUnit = string(unit);
                rock.perm.value = value * Convert(unit);           
                break
            elseif(strcmp(line(1),strcat('INITIALWATER')))
                rock.Swi.include = true;               
                rock.Swi.inputValue = value;
                rock.Swi.value = value;
                break
            elseif(strcmp(line(1),strcat('HETEROGENEOUS')))
                rock.heterogeneous = true;               
                break
            elseif(strcmp(line(1),strcat('POROSITY_PROFILE')))
                rock.poro.include = true;               
                rock.poro.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(rock.poro.fullFile);
                rock.poro.filePath = filePath;
                rock.poro.fileName = strcat(fileName,fileExt);
                rock.poro.porosity_profile = xlsread(rock.poro.fullFile);
                break
            elseif(strcmp(line(1),strcat('ALPHA')))
                rock.alpha.include = true;               
                rock.alpha.inputValue = value;
                rock.alpha.value = value;  
                break        
            elseif(strcmp(line(1),strcat('INDEX_MASK')))
                rock.het_index_mask.include = true;   
                C = line{2};
                rock.het_index_mask = str2double(regexp(C,'[\d.]+','match'));
                break    
            else
                break
            end
        end
    end
    experiment.rock = rock;
    
    %% block FLUID = ["DENSITYW" "DENSITYNW" "VISCOSITYW" "VISCOSITYNW"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            value = str2double(line(2));
            if(length(line) == 3), unit = line(3); end
            if(strcmp(line(1),strcat('DENSITYW')))
                fluid.rhoW.include = true;               
                fluid.rhoW.inputValue = value;
                fluid.rhoW.inputUnit = string(unit);
                fluid.rhoW.value = value * Convert(unit);         
                break         
            elseif(strcmp(line(1),strcat('DENSITYNW')))
                fluid.rhoNW.include = true;               
                fluid.rhoNW.inputValue = value;
                fluid.rhoNW.inputUnit = string(unit);
                fluid.rhoNW.value = value * Convert(unit);            
                break
            elseif(strcmp(line(1),strcat('VISCOSITYW')))
                fluid.muW.include = true;               
                fluid.muW.inputValue = value;
                fluid.muW.inputUnit = string(unit);
                fluid.muW.value = value * Convert(unit);       
                break
            elseif(strcmp(line(1),strcat('VISCOSITYNW')))
                fluid.muNW.include = true;               
                fluid.muNW.inputValue = value;
                fluid.muNW.inputUnit = string(unit);
                fluid.muNW.value = value * Convert(unit);              
                break
            else
                break
            end
        end
    end
    experiment.fluid = fluid;
    
    %% block PROCESS = ["SS" "USS" "CENT"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            type = line{1};
            name = line{2};
            if(strcmpi(line(1),strcat('SS')))
                process.type = type;
                process.name = name;                
                break
            elseif(strcmpi(line(1),strcat('USS')))
                process.type = type;
                process.name = name;
                break
            elseif(strcmpi(line(1),strcat('CENT')))
                process.type = type;
                process.name = name;
                break
            else
                break
            end
        end
    end
    process.include = true;
    experiment.process = process;
    
    %% block SIMULATION = ["TYPE" "NCELLS" "BCELLS" "TIMESTEP" "RAMPUPSTEPS" "GAUGEOFF"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    
    simulation.load_from_sat_prof = false;
    simulation.high_percision_mode = false; 

    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
         while (~isempty(find(~cellfun(@isempty,line),1)) && ...
                ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('TYPE')))
                simulation.type.include = true;
                simulation.type = line{2};
                break
            end    
            if length(line)>1; value = str2double(line(2)); end
            if(length(line) == 3), unit = line(3); end
            if(strcmp(line(1),strcat('NCELLS')))
                simulation.nCells.include = true;
                simulation.nCells.inputValue = value;
                simulation.nCells.value = value;
                break
            elseif(strcmp(line(1),strcat('BCELLS')))
                simulation.bCells.include = true;
                simulation.bCells.inputValue = value;
                simulation.bCells.firstCellSw.inputValue = str2double(line(3));
                simulation.bCells.lastCellSw.inputValue = str2double(line(4));   
                simulation.bCells.value = value;
                simulation.bCells.firstCellSw.value = str2double(line(3));
                simulation.bCells.lastCellSw.value = str2double(line(4));           
                break
            elseif(strcmp(line(1),strcat('MAXTIMESTEP')))
                simulation.timeStep.include = true;               
                simulation.timeStep.inputValue = value;
                simulation.timeStep.inputUnit = string(AppUnit(unit));
                simulation.timeStep.value = value * Convert(unit);      
                break   
            elseif(strcmp(line(1),strcat('RAMPUPSTEPS')))
                simulation.rampupsteps.include = true;               
                simulation.rampupsteps.inputValue = value;
                simulation.rampupsteps.value = value;      
                break 
            elseif(strcmp(line(1),strcat('LOAD_FROM_SAT_PROF')))
                simulation.load_from_sat_prof = true;                    
                break     
            elseif(strcmp(line(1),strcat('HIGH_PERCISION_MODE')))
                simulation.high_percision_mode = true;                    
                break 
            elseif(strcmp(line(1),strcat('GAUGEOFF')))
                simulation.gaugeOff.include = true;               
                simulation.gaugeOff.inputValue = value;
                simulation.gaugeOff.inputUnit = string(unit);
                simulation.gaugeOff.value = value * Convert(unit);
                break
            else
                break
            end      
        end
    end
    model.simulation = simulation;
        
    %% block SCHEDULE = ["FILENAME" "PINIT" "POUT" "CENTRAD" "STARTUP"]
    schedule = [];
    procedure = []; 
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('FILENAME')))
                schedule.include = true;
                schedule.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(schedule.fullFile);
                schedule.filePath = filePath;
                schedule.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(schedule.fullFile);
%                 [~,inputTable,~,table] = ImportTable2(schedule.fullFile);
%                 schedule.inputTable = inputTable;
                schedule.table = table;
                procedure = CreateProcedure(table, process.type);
                schedule.procedure = procedure;                
                break
            end       
            value = str2double(line(2));
            if(length(line) == 3), unit = line(3); end
            if(strcmp(line(1),strcat('PINI')))
                schedule.pini.include = true;               
                schedule.pini.inputValue = value;
                schedule.pini.inputUnit = string(unit);
                schedule.pini.value = value * Convert(unit);       
                break
            elseif(strcmp(line(1),strcat('POUT')))
                schedule.pout.include = true;               
                schedule.pout.inputValue = value;
                schedule.pout.inputUnit = string(unit);
                schedule.pout.value = value * Convert(unit);            
                break
            elseif(strcmp(line(1),strcat('CENTRAD'))) 
                if(strcmpi(process.type,'CENT'))
                    schedule.centRad.inputValue = value;
                    schedule.centRad.inputUnit = string(unit);
                    schedule.centRad.value = value * Convert(unit);
                    headers = procedure.Properties.VariableNames;
                    procedure_array = table2array(procedure);
                    procedure_array(:,3) = procedure_array(:,3) .^ 2 * schedule.centRad.value;
                    headers(3) = cellstr('Rotational Acceleration [m/s^2]');
                    procedure = array2table(procedure_array,'VariableNames',headers);
                    schedule.procedure = procedure;            
                end
                break
            elseif(strcmp(line(1),strcat('STARTUP')))
                if(strcmpi(process.type,'CENT'))
                    schedule.startup.include = true;                               
                    value = str2double(line(2));
                    unit = line(3);
                    schedule.startupPeriod.inputValue = value;
                    schedule.startupPeriod.inputUnit = string(unit);
                    schedule.startupPeriod.value = value * Convert(unit);
                    value_RPM = str2double(line(4));
                    unit_RPM = line(5);
                    schedule.startupRPM.inputValue = value_RPM;
                    schedule.startupRPM.inputUnit = string(unit_RPM);  
                    schedule.startupRPM.value = (value_RPM * Convert(unit_RPM)) ^ 2 * schedule.centRad.value;
                    headers = procedure.Properties.VariableNames;
                    procedure_array = table2array(procedure); 
                    schedule_in_rpm = sqrt(procedure_array(:,3)/schedule.centRad.value)/Convert('rpm');
                    procedure_table_in_rpm = [procedure_array(:,1:2),schedule_in_rpm];
                    start_up_rpm = sqrt(schedule.startupRPM.value/schedule.centRad.value)/Convert('rpm');
                    
                    startupSlope = start_up_rpm / ...
                                   schedule.startupPeriod.value;
                                        
                    added_schedule_row = [];
                    interval_between_steps = 20; % seconds
                    for i = 1 : height(procedure_table_in_rpm)
                        if i == 1
                            time_needed = procedure_table_in_rpm(i,3) / startupSlope;
                        else
                            time_needed = ( procedure_table_in_rpm(i,3) - procedure_table_in_rpm(i-1,3)) / startupSlope; 
                        end
                        residual_value = rem(time_needed, interval_between_steps);
                        number_of_steps = floor(time_needed / interval_between_steps);
                        added_schedule_times = linspace(procedure_table_in_rpm(i,1),procedure_table_in_rpm(i,1) + number_of_steps * interval_between_steps, number_of_steps + 1);
                        if i == 1
                            added_schedule_rpm = linspace(0 ,number_of_steps * interval_between_steps, number_of_steps + 1) * startupSlope;
                        else
                            added_schedule_rpm = linspace(0 ,number_of_steps * interval_between_steps, number_of_steps + 1) * startupSlope + procedure_table_in_rpm(i-1,3);   
                        end
                        for j = 1 : length(added_schedule_times) - 1
                            added_schedule_row = [added_schedule_row; ...
                                [added_schedule_times(j), added_schedule_times(j + 1), added_schedule_rpm(j + 1)]];
                        end
                        if not(residual_value == 0)
                            added_schedule_row = [added_schedule_row; ...
                                [added_schedule_times(end), added_schedule_times(end) + residual_value, added_schedule_rpm(end) + residual_value * startupSlope]];
                        end
                        added_schedule_row = [added_schedule_row; ...
                            [added_schedule_times(end) + residual_value, procedure_table_in_rpm(i,2), procedure_table_in_rpm(i,3)]];
                    end   


                    procedure = [added_schedule_row(:,1:2),(added_schedule_row(:,3)*Convert('rpm')) .^2 .* schedule.centRad.value];               
                    procedure = array2table(procedure,'VariableNames',headers);
                    schedule.procedure = procedure;
                end
                break                
            else
                break
            end                        
        end
    end
    experiment.schedule = schedule;
    experiment.schedule_input = schedule; %used for reversion
    
    %% block OBSERVATION = ["PRESSURE" "SWAVG" "SATPROFILE" "PRODUCTION"]
    observation = [];
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('PRESSURE')))
                observation.pressure.include = true;
                observation.pressure.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(observation.pressure.fullFile);
                observation.pressure.filePath = filePath;
                observation.pressure.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.pressure.fullFile);
%                 [~,inputTable,~,table] = ImportTable2(observation.pressure.fullFile);
%                 observation.pressure.inputTable = inputTable;
                observation.pressure.table = table;           
                break 
            elseif(strcmp(line(1),strcat('SWAVG')))
                observation.swavg.include = true;                
                observation.swavg.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(observation.swavg.fullFile);
                observation.swavg.filePath = filePath;
                observation.swavg.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.swavg.fullFile);
                observation.swavg.table = table;           
                break 
            elseif(strcmp(line(1),strcat('SATPROFILE')))
                observation.satProfile.include = true;
                observation.satProfile.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(observation.satProfile.fullFile);
                observation.satProfile.filePath = filePath;
                observation.satProfile.fileName = strcat(fileName,fileExt);                
                observation.satProfile.table = table2array(readtable(observation.satProfile.fullFile));
                observation.satProfile.table(1,1) = 0;             
                break 
            elseif(strcmp(line(1),strcat('PRODUCTION')))
                observation.prod.include = true;
                observation.prod.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(observation.prod.fullFile);
                observation.prod.filePath = filePath;
                observation.prod.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.prod.fullFile);
                observation.prod.table = table;
                break
            else
                break
            end
        end
    end
    experiment.observation = observation;  
    
    %% block SATURATION FUNCTIONS = ["PC" "KRW" "KRNW" "SWC" "SOR" "KRW@SOR" "KRO@SWC"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;    
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('KR')))                
                kr.type = string(line(2));
                kr.include = true;
                kr = read_kr_parameters_from_settings(kr, line);
                satfun.kr = kr;
                break
            end
            if(strcmp(line(1),strcat('PC')))
                pc.type = string(line{2});
                pc.include = true;
                pc = read_pc_parameters_from_settings(pc, line);
                satfun.pc = pc;
                break
            end
            if(strcmpi(line(1),strcat('KR_compare_1')))
                kr_compare_1.type = string(line(2));
                kr_compare_1.include = true;
                kr_compare_1 = read_kr_parameters_from_settings(kr_compare_1, line);
                satfun.kr_compare_1 = kr_compare_1;
                break
            end
            if(strcmpi(line(1),strcat('KR_ub')))
                kr_ub.type = string(line(2));
                kr_ub.include = true;
                kr_ub = read_kr_parameters_from_settings(kr_ub, line);
                satfun.kr_ub = kr_ub;
                break
            end
            if(strcmpi(line(1),strcat('KR_lb')))
                kr_lb.type = string(line(2));
                kr_lb.include = true;
                kr_lb = read_kr_parameters_from_settings(kr_lb, line);
                satfun.kr_lb = kr_lb;
                break
            end
            if(strcmpi(line(1),strcat('PC_compare_1')))
                pc_compare_1.type = string(line{2});
                pc_compare_1.include = true;
                pc_compare_1 = read_pc_parameters_from_settings(pc_compare_1, line);
                satfun.pc_compare_1 = pc_compare_1;
                break
            end
            if(strcmpi(line(1),strcat('PC_ub')))
                pc_ub.type = string(line{2});
                pc_ub.include = true;
                pc_ub = read_pc_parameters_from_settings(pc_ub, line);
                satfun.pc_ub = pc_ub;
                break
            end
            if(strcmpi(line(1),strcat('PC_lb')))
                pc_lb.type = string(line{2});
                pc_lb.include = true;
                pc_lb = read_pc_parameters_from_settings(pc_lb, line);
                satfun.pc_lb = pc_lb;
                break
            end
        end
    end
    experiment.satfun = satfun;    
    model.experiment = experiment;
    %% block PLOT OPTIONS = ["STYLE" "COLORMAP" "DISPLAYTIME"
    %                        "DISPLAYLENGTH" "DISPLAYVOLUME" 
    %                        "DISPLAYPRESS" "DISPLAYRATE"]
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4; 
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('STYLE')))
                plot.style.include = true;
                plot.style.inputStyle = line{2};
                break
            elseif(strcmp(line(1),strcat('COLORMAP')))
                plot.colormap.include = true;
                plot.colormap.inputColormap = line{2};
                break           
            elseif(strcmp(line(1),strcat('DISPLAYTIME')))
                plot.displayTime.include = true;
                plot.displayTime.inputUnit = AppUnit(line{2});
                break  
            elseif(strcmp(line(1),strcat('DISPLAYLENGTH')))
                plot.displayLength.include = true;
                plot.displayLength.inputUnit = AppUnit(line{2});
                break
            elseif(strcmp(line(1),strcat('DISPLAYVOLUME')))
                plot.displayVolume.include = true;
                plot.displayVolume.inputUnit = AppUnit(line{2});
                break
            elseif(strcmp(line(1),strcat('DISPLAYPRESS')))
                plot.displayPress.include = true;
                plot.displayPress.inputUnit = AppUnit(line{2});
                break  
            elseif(strcmp(line(1),strcat('DISPLAYRATE')))
                plot.displayRate.include = true;
                plot.displayRate.inputUnit = AppUnit(line{2});
                break
            else
                break
            end   
        end
    end
    model.plot = plot;
    
    %% block OUTPUT OPTIONS = ["PATH" "SATPROFILE" "SAVEVIDEO" "SAVECONFIG"
    %                          "QUANTITIES" "SWAVG" "INJ" "PROD" "DELTAP"]
    output = [];
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('PATH')))
                output.path.include = true;
                output.path.dir = line{2};          
                break
            elseif(strcmp(line(1),strcat('SATPROFILE')))
                output.include = true;
                output.satProfile.include = true;
                output.satProfile.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(output.satProfile.fullFile);
                output.satProfile.filePath = filePath;
                output.satProfile.fileName = strcat(fileName,fileExt);        
                break                    
            elseif(strcmp(line(1),strcat('SAVECONFIG')))
                output.include = true;
                output.saveConfig.include = true;
                output.saveConfig.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(output.saveConfig.fullFile);
                output.saveConfig.filePath = filePath;
                output.saveConfig.fileName = strcat(fileName,fileExt);
                break
            elseif(strcmp(line(1),strcat('QUANTITIES')))
                output.include = true;
                output.quantities.include = true;
                output.quantities.fullFile = line{2};
                [filePath,fileName,fileExt] = fileparts(output.quantities.fullFile);
                output.quantities.filePath = filePath;
                output.quantities.fileName = strcat(fileName,fileExt);              
                break         
            elseif(strcmp(line(1),strcat('TIME')))
                output.time.include = true;
                break
            elseif(strcmp(line(1),strcat('SWAVG')))
                output.swavg.include = true;
                break
            elseif(strcmp(line(1),strcat('INJ')))
                output.inj.include = true;
                break
            elseif(strcmp(line(1),strcat('PROD')))
                output.prod.include = true;
                break
            elseif(strcmp(line(1),strcat('DELTAP')))
                output.deltaP.include = true;
                break
            else
                break
            end          
        end        
    end
    model.output = output;
%% block GENERAL HISTORY MATCH CONFIGURATIONS = ["EXCEL_FILE_NAME","EXCEL_FILE_PATH",
    %"MULTIPOINT","MULTIPOINT_RANDOM_NUMBERS","KR","KR_MODEL","PC","PC_MODEL",
    %"PDIFF_WEIGHT","SWAVG_WEIGHT","PROD_WEIGHT","SAT_PROFILE_WEIGHT","PDIF_ERROR",
    %"SWAVG_ERROR","PROD_ERROR","SAT_PROFILE_ERROR"]
    history_match = [];
    history_match.multi_point = false;
    history_match.kr.status = false;
    history_match.pc.status = false;
    
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('EXCEL_FILE_NAME')))
                history_match.hm_template_name = line{2};       
                break
            elseif(strcmp(line(1),strcat('EXCEL_FILE_PATH')))
                if strcmpi(simulation.type, 'historymatch')
                    history_match.hm_template_path = line{2}; 
                    if(isfield(history_match,'hm_template_name'))
                        assert(isfile(fullfile(history_match.hm_template_path,...
                                history_match.hm_template_name)), 'Could not find the history match template!')
                            data = readtable(fullfile(history_match.hm_template_path,...
                                history_match.hm_template_name),'ReadVariableNames',...
                                true,'ReadRowNames',true,'HeaderLines',0);
                            history_match.RowNames = data.Properties.RowNames;
                            history_match.imported_excel_table = data;
                            history_match.ndim = length(data.Properties.RowNames);
                            history_match.x0 = data{:,1}(:)';
                            history_match.lb = data{:,2}(:)';
                            history_match.ub = data{:,3}(:)';
                            % Sw should be given in ascending order
                            history_match.kr.Sw_hm = data{:,4}(~isnan(data{:,4}))';
                            % Sw should be given in ascending order
                            history_match.pc.Sw_hm = data{:,5}(~isnan(data{:,5}))';
                            % setting MCMC specific parameters
                            history_match.init_std = data{:,6}(~isnan(data{:,6}))';
                    end
                end
                break 
            elseif(strcmp(line(1),strcat('MULTIPOINT')))
                history_match.multi_point = true;      
                break
            elseif(strcmp(line(1),strcat('MULTIPOINT_RANDOM_NUMBERS')))
                history_match.multipoint_rnd_num = str2double(line{2});  
                break
            elseif(strcmp(line(1),strcat('KR')))
                history_match.kr.status = true;     
                break
            elseif(strcmp(line(1),strcat('KR_MODEL')))
                history_match.kr.type = line{2};     
                break
            elseif(strcmp(line(1),strcat('PC')))
                history_match.pc.status = true;   
                break
            elseif(strcmp(line(1),strcat('PC_MODEL')))
                history_match.pc.type = line{2};    
                break
            elseif(strcmp(line(1),strcat('PDIFF_WEIGHT')))
                history_match.pdiff_weight = str2double(line{2});   
                break
            elseif(strcmp(line(1),strcat('SWAVG_WEIGHT')))
                history_match.swavg_weight = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('PROD_WEIGHT')))
                history_match.prod_weight = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('SAT_PROFILE_WEIGHT')))
                history_match.sat_profile_weight = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('PDIF_ERROR')))
                history_match.pdiff_error = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('SWAVG_ERROR')))
                history_match.swavg_error = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('PROD_ERROR')))
                history_match.prod_error = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('SAT_PROFILE_ERROR')))
                history_match.sat_profile_error = str2double(line{2});    
                break
            else
                break
            end        
        end
    end
    model.history_match = history_match;
    
%% block GRADIENT BASED HISTORY MATCH CONFIGURATIONS = ["USE_PARALLEL","OPTIMALITY_TOLERANCE",
    %"STEP_TOLERANCE","MAX_FUNCTION_EVALUATIONS","SCALE_PROBLEM","OBJECTIVE_FUNCTION_TYPE",
    %"CENT_FILE_NAME","CENT_FILE_PATH","HISTORYMATCH_ALGORITHM"]
    
    history_match.UseParallel = false;
    history_match.ScaleProblem = false;
    
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = blockStart(blockNo + 1) - 4;
%     lineEnd   = length(lines) - whiteSpace;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))        
            if(strcmp(line(1),strcat('USE_PARALLEL')))
                history_match.UseParallel = true;   
                break
            elseif(strcmp(line(1),strcat('OPTIMALITY_TOLERANCE')))
                history_match.OptimalityTolerance = str2double(line{2});  
                break
            elseif(strcmp(line(1),strcat('STEP_TOLERANCE')))
                history_match.StepTolerance = str2double(line{2});
                break
            elseif(strcmp(line(1),strcat('MAX_FUNCTION_EVALUATIONS')))
                history_match.MaxFunctionEvaluations = str2double(line{2});
                break
            elseif(strcmp(line(1),strcat('SCALE_PROBLEM')))
                history_match.ScaleProblem = true;   
                break
            elseif(strcmp(line(1),strcat('OBJECTIVE_FUNCTION_TYPE')))
                history_match.obj_fun = line{2}; 
                break
            elseif(strcmp(line(1),strcat('CENT_FILE_NAME')))
                history_match.Cent_file_name = line{2};    
                break
            elseif(strcmp(line(1),strcat('CENT_FILE_PATH')))
                history_match.Cent_file_path = line{2};  
                break
            elseif(strcmp(line(1),strcat('HISTORYMATCH_ALGORITHM')))
                history_match.algorithm = line{2};   
                break
            else
                break
            end
        end
    end
    model.history_match = history_match;
    
%% block MCMC SPECIFIC CONFIGURATIONS = ["RANDOM_SEED",
    %"SAMPLE_REFINEMENT_COUNT","CHAINSIZE","MPI_ENABLED"]
     
    history_match.mpi_enabled = false;
    
    blockNo = blockNo + 1;
    lineBegin = blockStart(blockNo);
    lineEnd   = length(lines) - whiteSpace;
    for k = lineBegin : lineEnd            
        line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
        while (~isempty(find(~cellfun(@isempty,line),1)) && ...
               ~startsWith(line{1},'#'))
            if(strcmp(line(1),strcat('RANDOM_SEED')))
                history_match.random_seed = str2double(line{2});    
                break
            elseif(strcmp(line(1),strcat('SAMPLE_REFINEMENT_COUNT')))
                history_match.sampleRefinementCount = str2double(line{2});   
                break
            elseif(strcmp(line(1),strcat('CHAINSIZE')))
                history_match.chainSize = str2double(line{2});
                break
            elseif(strcmp(line(1),strcat('MPI_ENABLED')))
                history_match.mpi_enabled = true;
                break
            else
                break
            end
        end
    end
    model.history_match = history_match;

    %% close the file
    fclose(file);  
end