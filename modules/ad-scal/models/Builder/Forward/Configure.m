function model = Configure(path_to_settings)
%
% DESCRIPTION: Reads the keywords and the corresponding values from the 
%              input settings file
%
% SYNOPSIS:
%   model = Configure(path, fileName)
%
% PARAMETERS:
%   path_to_file - string path to the directory of the settings file
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
[filepath,name] = fileparts(path_to_settings);

fprintf("Configuring from: %s\n", name);
fprintf("Located at: %s\n", filepath);

[fid, msg] = fopen(path_to_settings, 'rt');
if fid < 0, error([path_to_settings + ': ' + msg]), end

% initialization of structs
geometry = struct;
rock.heterogeneous = false;
fluid = struct;
process = struct;
schedule = struct;
experiment = struct;
simulation = struct;
output = [];
history_match = struct;
history_match.multi_point = false;
history_match.kr.status = false;
history_match.pc.status = false;
history_match.UseParallel = false;
history_match.ScaleProblem = false;

% going through the file
while ~feof(fid)
    lin = fgetl(fid);
    if lin == -1
        msg = ferror(fid, 'clear');
        fclose(fid);
        error('SCAL_Settings:Input:Empty', ...
            'Settings file ''%s'' is unreadable.\nSystem reports: %s\n', ...
            fn, msg);
    end
    % Loop until next keyword
    kw = regexp(lin, '^[A-Z][A-Z0-9_]{0,30}(|/)', 'match', 'once');
    while isempty(kw) && ~feof(fid)
      lin = fgetl(fid);
      if lin ~= -1
         kw = regexp(lin, '^[A-Z][A-Z0-9_]{0,30}(|/)', 'match', 'once');
      end
    end

    if ~feof(fid)
        splitted_lin = split(lin);
        switch kw
            %% block GEOMETRY = ["LENGTH" "DIAMETER"]    
            case 'LENGTH'
                geometry.length.include = true;
                value = str2double(splitted_lin{2});
                geometry.length.inputValue = value;
                unit = string(splitted_lin{3});
                geometry.length.inputUnit = unit;
                geometry.length.value = value * Convert(unit);
            case 'DIAMETER'
                geometry.diameter.include = true;
                value = str2double(splitted_lin{2});
                geometry.diameter.inputValue = value;
                unit = string(splitted_lin{3});
                geometry.diameter.inputUnit = string(unit);
                geometry.diameter.value = value * Convert(unit);
            %% block ROCK = ["POROSITY" "PERMEABILITY" "INITIALWATER"
            case 'POROSITY'
                value = str2double(splitted_lin{2});
                rock.poro.include = true;
                rock.poro.inputValue = value;
                rock.poro.value = value;
            case 'PERMEABILITY'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                rock.perm.include = true;
                rock.perm.inputValue = value;
                rock.perm.inputUnit = string(unit);
                rock.perm.value = value * Convert(unit); 
            case 'INITIALWATER'
                value = str2double(splitted_lin{2});
                rock.Swi.include = true;               
                rock.Swi.inputValue = value;
                rock.Swi.value = value;
            % Enables simulations including capillary heterogeneities/ !not
            % implemented yet
            case 'HETEROGENEOUS'
                rock.heterogeneous = true;
            % Read the porosity profile for heterogeneous modelling
            case 'POROSITY_PROFILE'
                rock.poro.include = true;
                rock.poro.fullFile = abs_path(splitted_lin{2});
                [filePath,fileName,fileExt] = fileparts(rock.poro.fullFile);
                rock.poro.filePath = filePath;
                rock.poro.fileName = strcat(fileName,fileExt);
                rock.poro.porosity_profile = readmatrix(rock.poro.fullFile);
            % Defines the trade off between pressure and saturation for
            % calculating f factors in heterogeneous modelling
            case 'ALPHA'
                value = str2double(splitted_lin{2});
                rock.alpha.include = true;               
                rock.alpha.inputValue = value;
                rock.alpha.value = value; 
            % Index of the saturation profile rows 
            case 'INDEX_MASK'
                rock.het_index_mask.include = true;   
                C = splitted_lin{2};
                rock.het_index_mask = str2double(regexp(C,'[\d.]+','match'));
            %% block FLUID = ["DENSITYW" "DENSITYNW" "VISCOSITYW" "VISCOSITYNW"]
            case 'DENSITYW'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                fluid.rhoW.include = true;               
                fluid.rhoW.inputValue = value;
                fluid.rhoW.inputUnit = string(unit);
                fluid.rhoW.value = value * Convert(unit);
            case 'DENSITYNW'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                fluid.rhoNW.include = true;               
                fluid.rhoNW.inputValue = value;
                fluid.rhoNW.inputUnit = string(unit);
                fluid.rhoNW.value = value * Convert(unit); 
            case 'VISCOSITYW'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                fluid.muW.include = true;               
                fluid.muW.inputValue = value;
                fluid.muW.inputUnit = string(unit);
                fluid.muW.value = value * Convert(unit); 
            case 'VISCOSITYNW'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                fluid.muNW.include = true;               
                fluid.muNW.inputValue = value;
                fluid.muNW.inputUnit = string(unit);
                fluid.muNW.value = value * Convert(unit); 
            %% block PROCESS = ["SS" "USS" "CENT"]
            case {'SS', 'USS', 'CENT'}
                type = string(splitted_lin{1});
                % name can be Drainage or Imbibition
                name = string(splitted_lin{2});
                process.type = type;
                process.name = name;
            %% block SIMULATION = ["TYPE" "NCELLS" "BCELLS" "MAXTIMESTEP"
            % "RAMPUPSTEPS" "LOAD_FROM_SAT_PROF" "HIGH_PERCISION_MODE" "GAUGEOFF"]
            % forward or history match
            case 'TYPE'
                simulation.type.include = true;
                simulation.type = string(splitted_lin{2});
            % number of grid cells
            case 'NCELLS'
                value = str2double(splitted_lin{2});
                simulation.nCells.include = true;
                simulation.nCells.inputValue = value;
                simulation.nCells.value = value;
            % number of boundary cells including their wetting phase saturations
            case 'NCELLS_Y'
                value = str2double(splitted_lin{2});
                simulation.nCells_y.include = true;
                simulation.nCells_y.inputValue = value;
                simulation.nCells_y.value = value;
            % boundary cells of the simulation with pc = 0
            case 'BCELLS'
                simulation.bCells.include = true;
                value = str2double(splitted_lin{2});
                simulation.bCells.inputValue = value;
                simulation.bCells.firstCellSw.inputValue = str2double(splitted_lin{3});
                simulation.bCells.lastCellSw.inputValue = str2double(splitted_lin{4});   
                simulation.bCells.value = value;
                simulation.bCells.firstCellSw.value = str2double(splitted_lin{3});
                simulation.bCells.lastCellSw.value = str2double(splitted_lin{4});
            % maximum time step used for MRST rampup function
            case 'MAXTIMESTEP'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                simulation.timeStep.include = true;               
                simulation.timeStep.inputValue = value;
                simulation.timeStep.inputUnit = string(AppUnit(unit));
                simulation.timeStep.value = value * Convert(unit); 
            % number of steps used in the MRST rampup function
            case 'RAMPUPSTEPS'
                value = str2double(splitted_lin{2});
                simulation.rampupsteps.include = true;               
                simulation.rampupsteps.inputValue = value;
                simulation.rampupsteps.value = value;  
            % if this one is activated the time steps are read from the
            % saturation profile and other time stepping settings are ignored
            case 'LOAD_FROM_SAT_PROF'
                simulation.load_from_sat_prof = true;  
            % sets a limit of 5% saturation change for the solver, above which
            % the time steps are cut to increase the accuracy of the results
            % but increases the simulation time 
            case 'HIGH_PERCISION_MODE'
               simulation.high_percision_mode = true;  
            % read the pressure not from the both ends but somewhere in the
            % middle defined using this keyword
            case 'GAUGEOFF'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                simulation.gaugeOff.include = true;               
                simulation.gaugeOff.inputValue = value;
                simulation.gaugeOff.inputUnit = string(unit);
                simulation.gaugeOff.value = value * Convert(unit);
            %% block SCHEDULE = ["FILENAME" "PINIT" "POUT" "CENTRAD" "STARTUP"]
            % reads the schedule file table
            case 'FILENAME'
                schedule.include = true;
                schedule.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(schedule.fullFile);
                schedule.filePath = filePath;
                schedule.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(schedule.fullFile);
                schedule.table = table;
                procedure = CreateProcedure(table, process.type);
                schedule.procedure = procedure; 
            % setup inlet initial pressure
            case 'PINI'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                schedule.pini.include = true;               
                schedule.pini.inputValue = value;
                schedule.pini.inputUnit = string(unit);
                schedule.pini.value = value * Convert(unit); 
            % outlet initial pressure
            case 'POUT'
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                schedule.pout.include = true;               
                schedule.pout.inputValue = value;
                schedule.pout.inputUnit = string(unit);
                schedule.pout.value = value * Convert(unit);
            % distance between rotation axis of the centrifuge and middle of
            % the core
            case 'CENTRAD'
                switch process.type
                    case 'CENT'
                        value = str2double(splitted_lin{2});
                        unit = string(splitted_lin{3});
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
            % start up time of the centrifuge experiment
            case 'STARTUP'
                switch process.type
                    case 'CENT'
                        value = str2double(splitted_lin{2});
                        unit = string(splitted_lin{3});
                        schedule.startup.include = true;
                        schedule.startupPeriod.inputValue = value;
                        schedule.startupPeriod.inputUnit = string(unit);
                        schedule.startupPeriod.value = value * Convert(unit);
                        value_RPM = str2double(splitted_lin{4});
                        unit_RPM = string(splitted_lin{5});
                        [procedure, schedule] = setup_startup(schedule,...
                            procedure, value_RPM, unit_RPM);
                        schedule.procedure = procedure;
                end
            %% block OBSERVATION = ["PRESSURE" "PRESSURE_MID" "PRESSURE_MID_GAUGEOFF" 
            % "SWAVG" "SATPROFILE" "PRODUCTION"]
            case 'PRESSURE'
                observation.pressure.include = true;
                observation.pressure.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(observation.pressure.fullFile);
                observation.pressure.filePath = filePath;
                observation.pressure.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.pressure.fullFile);
                observation.pressure.table = table;   
            % in case during an experiment we have pressure measurements from
            % both end and middle of the core, one can input them using this
            % keyword
            case 'PRESSURE_MID'
                observation.pressure_mid.include = true;
                observation.pressure_mid.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(observation.pressure_mid.fullFile);
                observation.pressure_mid.filePath = filePath;
                observation.pressure_mid.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.pressure_mid.fullFile);
                observation.pressure_mid.table = table; 
            % the place of the pressure tabs for the pressure measurements in
            % the middle of the core, in case of having two pressures at the
            % same time
            case 'PRESSURE_MID_GAUGEOFF'
                observation.pressure_mid.gaugeOff.include = true; 
                value = str2double(splitted_lin{2});
                unit = string(splitted_lin{3});
                observation.pressure_mid.gaugeOff.inputValue = value;
                observation.pressure_mid.gaugeOff.inputUnit = unit;
                observation.pressure_mid.gaugeOff.value = value * Convert(splitted_lin{3});
            case 'SWAVG'
                observation.swavg.include = true;                
                observation.swavg.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(observation.swavg.fullFile);
                observation.swavg.filePath = filePath;
                observation.swavg.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.swavg.fullFile);
                observation.swavg.table = table;   
            case 'SATPROFILE'
                observation.satProfile.include = true;
                observation.satProfile.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(observation.satProfile.fullFile);
                observation.satProfile.filePath = filePath;
                observation.satProfile.fileName = strcat(fileName,fileExt);                
                observation.satProfile.table = readmatrix(observation.satProfile.fullFile);
                observation.satProfile.table(1,1) = 0;
            case 'PRODUCTION'
                observation.prod.include = true;
                observation.prod.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(observation.prod.fullFile);
                observation.prod.filePath = filePath;
                observation.prod.fileName = strcat(fileName,fileExt);
                [~,table] = ImportTable(observation.prod.fullFile);
                observation.prod.table = table;
            %% block SATURATION FUNCTIONS = ["PC" "KR" "KR_compare" "KR_BOUNDARY"
            % "PC_COMPARE" "PC_BOUNDARY" "CLOSURE_CORR"]
            case 'KR'
                kr.type = string(splitted_lin{2});
                kr.include = true;
                kr = read_kr_parameters_from_settings(kr, splitted_lin, path_to_settings);
                satfun.kr = kr;
            case 'PC'
                pc.type = string(splitted_lin{2});
                pc.include = true;
                pc = read_pc_parameters_from_settings(pc, splitted_lin, path_to_settings);
                satfun.pc = pc;
            % in case we want to compare the input kr with another one, the
            % keyword should be used as "KR_COMPARE_1" or "KR_COMPARE_2" ...
            % to be able to input up to 4 curves
            case {'KR_COMPARE_1', 'KR_COMPARE_2', 'KR_COMPARE_3', 'KR_COMPARE_4'}
                l1_string = char(splitted_lin{1}); 
                index = str2double(regexp(l1_string,'\d*','match'));
                kr_compare{index}.type = string(splitted_lin{2});
                kr_compare{index}.include = true;
                kr_compare{index} = read_kr_parameters_from_settings(kr_compare{index}, splitted_lin);
                satfun.kr_compare{index} = kr_compare{index};
            % boundary to compare with the model input
            % keyword should be used as "KR_BOUNDARY_1" or "KR_BOUNDARY_2" ...
            case {'KR_BOUNDARY_1', 'KR_BOUNDARY_2', 'KR_BOUNDARY_3', 'KR_BOUNDARY_4'}
                l1_string = char(splitted_lin{1}); 
                index = str2double(regexp(l1_string,'\d*','match'));
                kr_boundary{index}.type = string(splitted_lin{2});
                kr_boundary{index}.include = true;
                kr_boundary{index} = read_kr_parameters_from_settings(kr_boundary{index}, splitted_lin);
                satfun.kr_boundary{index} = kr_boundary{index};
            % same as "KR_COMPARE" but for pc
            case {'PC_COMPARE_1', 'PC_COMPARE_2', 'PC_COMPARE_3', 'PC_COMPARE_4'}
                l1_string = char(splitted_lin{1}); 
                index = str2double(regexp(l1_string,'\d*','match'));
                pc_compare{index}.type = string(splitted_lin{2});
                pc_compare{index}.include = true;
                pc_compare{index} = read_pc_parameters_from_settings(pc_compare{index}, splitted_lin);
                satfun.pc_compare{index} = pc_compare{index};
            % same as "KR_BOUNDARY" but for pc
            case {'PC_BOUNDARY_1', 'PC_BOUNDARY_2', 'PC_BOUNDARY_3', 'PC_BOUNDARY_4'}
                l1_string = char(splitted_lin{1}); 
                pc_boundary{index}.type = string(splitted_lin{2});
                pc_boundary{index}.include = true;
                pc_boundary{index} = read_pc_parameters_from_settings(pc_boundary{index}, splitted_lin);
                satfun.pc_boundary{index} = pc_boundary{index};
            % applies a closure correction to the input pc, meaning the input
            % is in saturation units and the pc shifted towards water
            % saturation of 1 within this amount
            case 'CLOSURE_CORR'
                pc_closure.include = true;
                pc_closure.inputValue = str2double(splitted_lin{2});
                pc_closure.value = str2double(splitted_lin{2});
                satfun.pc_closure = pc_closure;
            %% block PLOT OPTIONS = ["STYLE" "COLORMAP" "DISPLAYTIME"
            %                        "DISPLAYLENGTH" "DISPLAYVOLUME" 
            %                        "DISPLAYPRESS" "DISPLAYRATE"]
            case 'STYLE'
                plot.style.include = true;
                plot.style.inputStyle = string(splitted_lin{2});
            case 'COLORMAP'
                plot.colormap.include = true;
                plot.colormap.inputColormap = string(splitted_lin{2});
            case 'DISPLAYTIME'
                plot.displayTime.include = true;
                plot.displayTime.inputUnit = AppUnit(string(splitted_lin{2}));
            case 'DISPLAYLENGTH'
                plot.displayLength.include = true;
                plot.displayLength.inputUnit = AppUnit(string(splitted_lin{2}));
            case 'DISPLAYVOLUME'
                plot.displayVolume.include = true;
                plot.displayVolume.inputUnit = AppUnit(string(splitted_lin{2}));
            case 'DISPLAYPRESS'
                plot.displayPress.include = true;
                plot.displayPress.inputUnit = AppUnit(string(splitted_lin{2}));
            case 'DISPLAYRATE'
                plot.displayRate.include = true;
                plot.displayRate.inputUnit = AppUnit(string(splitted_lin{2}));
            %% block OUTPUT OPTIONS = ["PATH" "SATPROFILE" "SAVEVIDEO" "SAVECONFIG"
            % "QUANTITIES" "SWAVG" "INJ" "PROD" "DELTAP"]
            % path to output saturation profile
            case 'SATPROFILE_OUT'
                output.include = true;
                output.satProfile.include = true;
                output.satProfile.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(output.satProfile.fullFile);
                output.satProfile.filePath = filePath;
                output.satProfile.fileName = strcat(fileName,fileExt); 
            % path to save the configuration file
            case 'SAVECONFIG'
                output.include = true;
                output.saveConfig.include = true;
                output.saveConfig.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(output.saveConfig.fullFile);
                output.saveConfig.filePath = filePath;
                output.saveConfig.fileName = strcat(fileName,fileExt);
            % indicate path to save output quantities
            case 'QUANTITIES'
                output.include = true;
                output.quantities.include = true;
                output.quantities.fullFile = abs_path(splitted_lin{2}, path_to_settings);
                [filePath,fileName,fileExt] = fileparts(output.quantities.fullFile);
                output.quantities.filePath = filePath;
                output.quantities.fileName = strcat(fileName,fileExt);
            % active to print in the output file
            case 'TIME'
                output.time.include = true;
            % active to print in the output file
            case 'SWAVG_OUT'
                output.swavg.include = true;
            % active to print in the output file
            case 'INJ'
                output.inj.include = true;
            % active to print in the output file
            case 'PROD'
                output.prod.include = true;
            % active to print in the output file
            case 'DELTAP'
                output.deltaP.include = true;
            %% block GENERAL HISTORY MATCH CONFIGURATIONS = ["EXCEL_FILE_NAME","EXCEL_FILE_PATH",
            %"MULTIPOINT","MULTIPOINT_RANDOM_NUMBERS","KR","KR_MODEL","PC","PC_MODEL",
            %"PDIFF_WEIGHT","SWAVG_WEIGHT","PROD_WEIGHT","SAT_PROFILE_WEIGHT","PDIF_ERROR",
            %"SWAVG_ERROR","PROD_ERROR","SAT_PROFILE_ERROR"]
            % name of the excel template used for inputing parameters for
            % the history matching  
            case 'EXCEL_FILE_NAME'
                history_match.hm_template_name = string(splitted_lin{2});
            % path of the excel file mentioned above
            case 'EXCEL_FILE_PATH'
                switch simulation.type
                    case 'historymatch'
                        history_match.hm_template_path = abs_path(splitted_lin{2}, path_to_settings);
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
            % activate multi point start history matching
            case 'MULTIPOINT'
                history_match.multi_point = true;
            % number of random starting points to start from
            case 'MULTIPOINT_RANDOM_NUMBERS'
                history_match.multipoint_rnd_num = str2double(splitted_lin{2}); 
            % activate kr history matching
            case 'KR_HM'
                history_match.kr.status = true;
            % define the model used for kr history matching
            case 'KR_MODEL'
                history_match.kr.type = string(splitted_lin{2});
            % activate pc history matching
            case 'PC_HM'
                history_match.pc.status = true;
            % define the model used for pc history matching
            case 'PC_MODEL'
                 history_match.pc.type = string(splitted_lin{2});
             % weight of the pressure error in the objective function
            case 'PDIFF_WEIGHT'
                history_match.pdiff_weight = str2double(splitted_lin{2});
            % weight of average water saturation error in the objective
            % function
            case 'SWAVG_WEIGHT'
                history_match.swavg_weight = str2double(splitted_lin{2});
            % weight of production error in the objective function
            case 'PROD_WEIGHT'
                history_match.prod_weight = str2double(splitted_lin{2});
            % weight of saturation profile error in the objective function
            case 'SAT_PROFILE_WEIGHT'
                history_match.sat_profile_weight = str2double(splitted_lin{2}); 
            % based on the objective functions defined, we use error 0 for the
            % gradient based history matching and an error percentage for each
            % measurement for the MCMC simulations
            case 'PDIF_ERROR'
                history_match.pdiff_error = str2double(splitted_lin{2});
            case 'SWAVG_ERROR'
                history_match.swavg_error = str2double(splitted_lin{2});
            case 'PROD_ERROR'
                history_match.prod_error = str2double(splitted_lin{2});
            case 'SAT_PROFILE_ERROR'
                history_match.sat_profile_error = str2double(splitted_lin{2});
            %% block GRADIENT BASED HISTORY MATCH CONFIGURATIONS = ["USE_PARALLEL","OPTIMALITY_TOLERANCE",
            %"STEP_TOLERANCE","MAX_FUNCTION_EVALUATIONS","SCALE_PROBLEM","OBJECTIVE_FUNCTION_TYPE",
            %"CENT_FILE_NAME","CENT_FILE_PATH","HISTORYMATCH_ALGORITHM"]
            case 'USE_PARALLEL'
                history_match.UseParallel = true;
            case 'OPTIMALITY_TOLERANCE'
                history_match.OptimalityTolerance = str2double(splitted_lin{2});  
            case 'STEP_TOLERANCE'
                history_match.StepTolerance = str2double(splitted_lin{2});
            case 'MAX_FUNCTION_EVALUATIONS'
                history_match.MaxFunctionEvaluations = str2double(splitted_lin{2});
            case 'SCALE_PROBLEM'
                history_match.ScaleProblem = true;   
            % use normal or simultaneous history matching here
            case 'OBJECTIVE_FUNCTION_TYPE'
                history_match.obj_fun = string(splitted_lin{2}); 
            % name of the centrifuge settings file for simultaneous history
            % matching
            case 'CENT_FILE_NAME'
                history_match.Cent_file_name = string(splitted_lin{2}); 
            % path of centrifuge settings file
            case 'CENT_FILE_PATH'
                history_match.Cent_file_path = abs_path(splitted_lin{2}, path_to_settings);  
            % can be all the algorithms in the fmincon function or
            % ga_multi_objective for multiobjective optimization using gene
            % algorithm
            case 'HISTORYMATCH_ALGORITHM'
                history_match.algorithm = string(splitted_lin{2});  
            %% block MCMC SPECIFIC CONFIGURATIONS = ["RANDOM_SEED",
            %"SAMPLE_REFINEMENT_COUNT","CHAINSIZE","MPI_ENABLED"]
            % random seed for MCMC simulation
            case 'RANDOM_SEED'
                history_match.random_seed = str2double(splitted_lin{2});                
            % refer to paramonte docs for this
            case 'SAMPLE_REFINEMENT_COUNT'
                history_match.sampleRefinementCount = str2double(splitted_lin{2});
            % number of required accpected samples in the MCMC chain
            case 'CHAINSIZE'
                history_match.chainSize = str2double(splitted_lin{2});
            % active for multi cpu computation, !use only when called using
            % mpiexe in the command prompt
            case 'MPI_ENABLED'
                history_match.mpi_enabled = true;
        end
    end
end
%% store the data in structures
experiment.geometry = geometry;
experiment.rock = rock;
experiment.fluid = fluid;
experiment.process = process;
experiment.schedule = schedule;
% used for reversion back to the input schedule without start-up 
experiment.schedule_input = schedule; 
experiment.observation = observation; 
experiment.satfun = satfun;  
model.plot = plot;
model.output = output;
model.experiment = experiment;
model.simulation = simulation; 
model.history_match = history_match;

%% close the file

fclose(fid);
% end of main function
end
%% Translate relative pathname to absolute pathname.
function path = abs_path(inc_fn, setting_path)
inc_fn(or(inc_fn == '/', inc_fn == '\')) = filesep;
    if inc_fn(1) ~= filesep
       % Translate relative pathname to absolute pathname.
       path = fullfile(fileparts(setting_path), inc_fn);
    end
end