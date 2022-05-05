function model = Configure(path, fileName)
%
% DESCRIPTION: Reads the keywords and the corresponding values from the 
%              input settings file
%
% SYNOPSIS:
%   model = Configure(path, fileName)
%
% PARAMETERS:
%   path - string path to the directory of the settings file
%   fileName - name of the settings file
%
% RETURNS:
%   model - struct containing the following fields:
%   - geometry: length and diameter of the core
%   - rock: rock properties like porosity and absolute permeability
%   - fluid: fluid properties like viscosities and densities
%   - process: type of the simulation - SS, USS or Centrifuge, drainage or
%   imbibition
%   - simulation: time stepping and grid cells information
%   - schedule: flooding or rotation schedule depending on the type of the
%   experiment
%   - observation: path to files containing the measured experimental data
%   used for comparison in plots or history matching
%   - experiment: saturation functions used for forward modeling
%   - plot: plotting options
%   - output: define the properties desired in the output
%   - history_match: gradient based and MCMC history matching options
%
% EXAMPLE:
%   model = Configure(".\","settings_case1.txt");
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
  
%% block ROCK = ["POROSITY" "PERMEABILITY" "INITIALWATER" 
% "HETEROGENEOUS" "POROSITY_PROFILE" "ALPHA" "INDEX_MASK"]
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
        % Enables simulations including capillary heterogeneities
        elseif(strcmp(line(1),strcat('HETEROGENEOUS')))
            rock.heterogeneous = true;               
            break
        % Read the porosity profile for heterogeneous modelling
        elseif(strcmp(line(1),strcat('POROSITY_PROFILE')))
            rock.poro.include = true;               
            rock.poro.fullFile = line{2};
            [filePath,fileName,fileExt] = fileparts(rock.poro.fullFile);
            rock.poro.filePath = filePath;
            rock.poro.fileName = strcat(fileName,fileExt);
            rock.poro.porosity_profile = xlsread(rock.poro.fullFile);
            break
        % Defines the trade off between pressure and saturation for
        % calculating f factors in heterogeneous modelling
        elseif(strcmp(line(1),strcat('ALPHA')))
            rock.alpha.include = true;               
            rock.alpha.inputValue = value;
            rock.alpha.value = value;  
            break     
        % Index of the saturation profile rows 
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
        % name can be Drainage or Imbibition
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

%% block SIMULATION = ["TYPE" "NCELLS" "BCELLS" "MAXTIMESTEP"
% "RAMPUPSTEPS" "LOAD_FROM_SAT_PROF" "HIGH_PERCISION_MODE" "GAUGEOFF"]
blockNo = blockNo + 1;
lineBegin = blockStart(blockNo);
lineEnd   = blockStart(blockNo + 1) - 4;

simulation.load_from_sat_prof = false;
simulation.high_percision_mode = false; 

for k = lineBegin : lineEnd            
    line = strtrim(regexp(lines{k},'\s+','split')); % split by spaces
     while (~isempty(find(~cellfun(@isempty,line),1)) && ...
            ~startsWith(line{1},'#'))
        % forward or history match
        if(strcmp(line(1),strcat('TYPE')))
            simulation.type.include = true;
            simulation.type = line{2};
            break
        end    
        % number of grid cells
        if length(line)>1; value = str2double(line(2)); end
        if(length(line) == 3), unit = line(3); end
        if(strcmp(line(1),strcat('NCELLS')))
            simulation.nCells.include = true;
            simulation.nCells.inputValue = value;
            simulation.nCells.value = value;
            break
        % number of boundary cells including their wetting phase saturations
        elseif(strcmp(line(1),strcat('NCELLS_Y')))
            simulation.nCells_y.include = true;
            simulation.nCells_y.inputValue = value;
            simulation.nCells_y.value = value;
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
        % maximum time step used for MRST rampup function
        elseif(strcmp(line(1),strcat('MAXTIMESTEP')))
            simulation.timeStep.include = true;               
            simulation.timeStep.inputValue = value;
            simulation.timeStep.inputUnit = string(AppUnit(unit));
            simulation.timeStep.value = value * Convert(unit);      
            break   
        % number of steps used in the MRST rampup function
        elseif(strcmp(line(1),strcat('RAMPUPSTEPS')))
            simulation.rampupsteps.include = true;               
            simulation.rampupsteps.inputValue = value;
            simulation.rampupsteps.value = value;      
            break 
        % if this one is activated the time steps are read from the
        % saturation profile and other time stepping settings are ignored
        elseif(strcmp(line(1),strcat('LOAD_FROM_SAT_PROF')))
            simulation.load_from_sat_prof = true;                    
            break     
        % sets a limit of 5% saturation change for the solver, above which
        % the time steps are cut to increase the accuracy of the results
        % but increases the simulation time
        elseif(strcmp(line(1),strcat('HIGH_PERCISION_MODE')))
            simulation.high_percision_mode = true;                    
            break 
        % read the pressure not from the both ends but somewhere in the
        % middle defined using this keyword
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
            schedule.table = table;
            procedure = CreateProcedure(table, process.type);
            schedule.procedure = procedure;                
            break
        end       
        % setup inlet initial pressure
        value = str2double(line(2));
        if(length(line) == 3), unit = line(3); end
        if(strcmp(line(1),strcat('PINI')))
            schedule.pini.include = true;               
            schedule.pini.inputValue = value;
            schedule.pini.inputUnit = string(unit);
            schedule.pini.value = value * Convert(unit);       
            break
        % outlet initial pressure
        elseif(strcmp(line(1),strcat('POUT')))
            schedule.pout.include = true;               
            schedule.pout.inputValue = value;
            schedule.pout.inputUnit = string(unit);
            schedule.pout.value = value * Convert(unit);            
            break
        % distance between rotation axis of the centrifuge and middle of
        % the core
        elseif(strcmp(line(1),strcat('CENTRAD'))) 
            if(strcmpi(process.type,'CENT'))
                schedule.centRad.inputValue = value;
                schedule.centRad.inputUnit = string(unit);
                schedule.centRad.value = value * Convert(unit);
                headers = procedure.Properties.VariableNames;
                procedure_array = table2array(procedure);
                procedure_array(:,3) = procedure_array(:,3) .^ 2 * ...
                    schedule.centRad.value;
                headers(3) = cellstr('Rotational Acceleration [m/s^2]');
                procedure = array2table(procedure_array,'VariableNames',...
                    headers);
                schedule.procedure = procedure;            
            end
            break
        % start up time of the centrifuge experiment
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
                schedule.startupRPM.value = (value_RPM * Convert(unit_RPM))...
                    ^ 2 * schedule.centRad.value;
                headers = procedure.Properties.VariableNames;
                procedure_array = table2array(procedure); 
                % schedule_in_rpm
                SIR = sqrt(procedure_array(:,3)/schedule.centRad.value)...
                    / Convert('rpm');
                % procedure_table_in_rpm
                PTIR = [procedure_array(:,1:2),SIR];
                % start_up_rpm
                SUR = sqrt(schedule.startupRPM.value/schedule.centRad.value)...
                    / Convert('rpm');
                % startupSlope
                stSL = SUR / schedule.startupPeriod.value;
                % added_schedule_row                 
                ASRow = [];
                % interval_between_steps
                IBS = 20; % seconds
                for i = 1 : height(PTIR)
                    if i == 1
                        % time_needed
                        TN = PTIR(i,3) / stSL;
                    else
                        TN = ( PTIR(i,3) - PTIR(i-1,3)) / stSL; 
                    end
                    % residual_value
                    RV = rem(TN, IBS);
                    % number_of_steps
                    NOS = floor(TN / IBS);
                    % added_schedule_times
                    AST = linspace(PTIR(i,1),PTIR(i,1) + NOS * IBS, NOS + 1);
                    if i == 1
                        % added_schedule_rpm
                        ASrpm = linspace(0 ,NOS * IBS, NOS + 1) * stSL;
                    else
                        ASrpm = linspace(0 ,NOS * IBS, NOS + 1) * stSL + ...
                            PTIR(i-1,3);   
                    end
                    for j = 1 : length(AST) - 1
                        ASRow = [ASRow; ...
                            [AST(j), AST(j + 1), ASrpm(j + 1)]];
                    end
                    if not(RV == 0)
                        ASRow = [ASRow; ...
                            [AST(end), AST(end) + RV, ASrpm(end) + RV * stSL]];
                    end
                    ASRow = [ASRow; ...
                        [AST(end) + RV, PTIR(i,2), PTIR(i,3)]];
                end   
                procedure = [ASRow(:,1:2),(ASRow(:,3)*Convert('rpm')) .^2 .* ...
                    schedule.centRad.value];               
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
% used for reversion back to the input schedule without start-up 
experiment.schedule_input = schedule; 

%% block OBSERVATION = ["PRESSURE" "PRESSURE_MID" "PRESSURE_MID_GAUGEOFF" 
% "SWAVG" "SATPROFILE" "PRODUCTION"]
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
            observation.pressure.table = table;           
            break 
        % in case during an experiment we have pressure measurements from
        % both end and middle of the core, one can input them using this
        % keyword
        elseif(strcmp(line(1),strcat('PRESSURE_MID')))
            observation.pressure_mid.include = true;
            observation.pressure_mid.fullFile = line{2};
            [filePath,fileName,fileExt] = fileparts(observation.pressure_mid.fullFile);
            observation.pressure_mid.filePath = filePath;
            observation.pressure_mid.fileName = strcat(fileName,fileExt);
            [~,table] = ImportTable(observation.pressure_mid.fullFile);
            observation.pressure_mid.table = table;           
            break 
        % the place of the pressure tabs for the pressure measurements in
        % the middle of the core, in case of having two pressures at the
        % same time
        elseif(strcmp(line(1),strcat('PRESSURE_MID_GAUGEOFF')))
            observation.pressure_mid.gaugeOff.include = true;               
            observation.pressure_mid.gaugeOff.inputValue = str2double(line{2});
            observation.pressure_mid.gaugeOff.inputUnit = string(line{3});
            observation.pressure_mid.gaugeOff.value = str2double(line{2}) * Convert(line{3});
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

%% block SATURATION FUNCTIONS = ["PC" "KR" "KR_compare" "KR_BOUNDARY"
% "PC_COMPARE" "PC_BOUNDARY" "CLOSURE_CORR"]
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
        elseif(strcmp(line(1),strcat('PC')))
            pc.type = string(line{2});
            pc.include = true;
            pc = read_pc_parameters_from_settings(pc, line);
            satfun.pc = pc;
            break
        % in case we want to compare the input kr with another one, the
        % keyword should be used as "KR_COMPARE_1" or "KR_COMPARE_2" ...
        % to be able to input as many as required
        elseif contains(line(1),"KR_COMPARE",'IgnoreCase',true)
            l1_string = char(line(1)); index = str2double(regexp(l1_string,'\d*','match'));
            kr_compare{index}.type = string(line(2));
            kr_compare{index}.include = true;
            kr_compare{index} = read_kr_parameters_from_settings(kr_compare{index}, line);
            satfun.kr_compare{index} = kr_compare{index};
            break
        % boundary to compare with the model input
        % keyword should be used as "KR_BOUNDARY_1" or "KR_BOUNDARY_2" ...
        elseif contains(line(1),"KR_BOUNDARY",'IgnoreCase',true)
            l1_string = char(line(1)); index = str2double(regexp(l1_string,'\d*','match'));
            kr_boundary{index}.type = string(line(2));
            kr_boundary{index}.include = true;
            kr_boundary{index} = read_kr_parameters_from_settings(kr_boundary{index}, line);
            satfun.kr_boundary{index} = kr_boundary{index};
            break
        % same as "KR_COMPARE" but for pc
        elseif contains(line(1),"PC_COMPARE",'IgnoreCase',true)
            l1_string = char(line(1)); index = str2double(regexp(l1_string,'\d*','match'));
            pc_compare{index}.type = string(line{2});
            pc_compare{index}.include = true;
            pc_compare{index} = read_pc_parameters_from_settings(pc_compare{index}, line);
            satfun.pc_compare{index} = pc_compare{index};
            break
        % same as "KR_BOUNDARY" but for pc
        elseif contains(line(1),"PC_BOUNDARY",'IgnoreCase',true)
            l1_string = char(line(1)); index = str2double(regexp(l1_string,'\d*','match'));
            pc_boundary{index}.type = string(line{2});
            pc_boundary{index}.include = true;
            pc_boundary{index} = read_pc_parameters_from_settings(pc_boundary{index}, line);
            satfun.pc_boundary{index} = pc_boundary{index};
            break
        % applies a closure correction to the input pc, meaning the input
        % is in saturation units and the pc shifted towards water
        % saturation of 1 within this amount
        elseif(strcmpi(line(1),strcat('CLOSURE_CORR')))
            pc_closure.include = true;
            pc_closure.inputValue = str2double(line(2));
            pc_closure.value = str2double(line(2));
            satfun.pc_closure = pc_closure;
            break
        else
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
        % path to save the requested output
        if(strcmp(line(1),strcat('PATH')))
            output.path.include = true;
            output.path.dir = line{2};          
            break
        % path to output saturation profile
        elseif(strcmp(line(1),strcat('SATPROFILE')))
            output.include = true;
            output.satProfile.include = true;
            output.satProfile.fullFile = line{2};
            [filePath,fileName,fileExt] = fileparts(output.satProfile.fullFile);
            output.satProfile.filePath = filePath;
            output.satProfile.fileName = strcat(fileName,fileExt);        
            break      
        % path to save the configuration file
        elseif(strcmp(line(1),strcat('SAVECONFIG')))
            output.include = true;
            output.saveConfig.include = true;
            output.saveConfig.fullFile = line{2};
            [filePath,fileName,fileExt] = fileparts(output.saveConfig.fullFile);
            output.saveConfig.filePath = filePath;
            output.saveConfig.fileName = strcat(fileName,fileExt);
            break
        % use alone to activate output quantities 
        elseif(strcmp(line(1),strcat('QUANTITIES')))
            output.include = true;
            output.quantities.include = true;
            output.quantities.fullFile = line{2};
            [filePath,fileName,fileExt] = fileparts(output.quantities.fullFile);
            output.quantities.filePath = filePath;
            output.quantities.fileName = strcat(fileName,fileExt);              
            break    
        % active to print in the output file
        elseif(strcmp(line(1),strcat('TIME')))
            output.time.include = true;
            break
        % active to print in the output file
        elseif(strcmp(line(1),strcat('SWAVG')))
            output.swavg.include = true;
            break
        % active to print in the output file
        elseif(strcmp(line(1),strcat('INJ')))
            output.inj.include = true;
            break
        % active to print in the output file
        elseif(strcmp(line(1),strcat('PROD')))
            output.prod.include = true;
            break
        % active to print in the output file
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
        % name of the excel template used for inputing parameters for
        % the history matching
        if(strcmp(line(1),strcat('EXCEL_FILE_NAME')))
            history_match.hm_template_name = line{2};       
            break
        % path of the excel file mentioned above
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
        % activate multi point start history matching
        elseif(strcmp(line(1),strcat('MULTIPOINT')))
            history_match.multi_point = true;      
            break
        % number of random starting points to start from
        elseif(strcmp(line(1),strcat('MULTIPOINT_RANDOM_NUMBERS')))
            history_match.multipoint_rnd_num = str2double(line{2});  
            break
        % activate kr history matching
        elseif(strcmp(line(1),strcat('KR')))
            history_match.kr.status = true;     
            break
       % define the model used for kr history matching
        elseif(strcmp(line(1),strcat('KR_MODEL')))
            history_match.kr.type = line{2};     
            break
        % activate pc history matching
        elseif(strcmp(line(1),strcat('PC')))
            history_match.pc.status = true;   
            break
        % define the model used for pc history matching
        elseif(strcmp(line(1),strcat('PC_MODEL')))
            history_match.pc.type = line{2};    
            break
        % weight of the pressure error in the objective function
        elseif(strcmp(line(1),strcat('PDIFF_WEIGHT')))
            history_match.pdiff_weight = str2double(line{2});   
            break
        % weight of average water saturation error in the objective
        % function
        elseif(strcmp(line(1),strcat('SWAVG_WEIGHT')))
            history_match.swavg_weight = str2double(line{2});    
            break
        % weight of production error in the objective function
        elseif(strcmp(line(1),strcat('PROD_WEIGHT')))
            history_match.prod_weight = str2double(line{2});    
            break
        % weight of saturation profile error in the objective function
        elseif(strcmp(line(1),strcat('SAT_PROFILE_WEIGHT')))
            history_match.sat_profile_weight = str2double(line{2});    
            break
        % based on the objective functions defined, we use error 0 for the
        % gradient based history matching and an error percentage for each
        % measurement for the MCMC simulations
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
        % active parallel cpu computation
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
        % use normal or simultaneous history matching here
        elseif(strcmp(line(1),strcat('OBJECTIVE_FUNCTION_TYPE')))
            history_match.obj_fun = line{2}; 
            break
        % name of the centrifuge settings file for simultaneous history
        % matching
        elseif(strcmp(line(1),strcat('CENT_FILE_NAME')))
            history_match.Cent_file_name = line{2};    
            break
        % path of centrifuge settings file
        elseif(strcmp(line(1),strcat('CENT_FILE_PATH')))
            history_match.Cent_file_path = line{2};  
            break
        % can be all the algorithms in the fmincon function or
        % ga_multi_objective for multiobjective optimization using genetic
        % algorithm
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
        % random seed for MCMC simulation
        if(strcmp(line(1),strcat('RANDOM_SEED')))
            history_match.random_seed = str2double(line{2});    
            break
        % refer to paramonte docs for this
        elseif(strcmp(line(1),strcat('SAMPLE_REFINEMENT_COUNT')))
            history_match.sampleRefinementCount = str2double(line{2});   
            break
        % number of required accpected samples in the MCMC chain
        elseif(strcmp(line(1),strcat('CHAINSIZE')))
            history_match.chainSize = str2double(line{2});
            break
        % active for multi cpu computation, !use only when called using
        % mpiexe in the command prompt
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