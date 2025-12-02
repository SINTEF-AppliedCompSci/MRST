function SaveConfig(model, fullFile)
%
% DESCRIPTION: saves the information in the model struct as a settings
% files which can be used in the Configure module - basically a way to save
% your work and start again
%
% SYNOPSIS:
%   SaveConfig(model, fullFile)
%
% PARAMETERS:
%   - model - main struct of the simulation
%   - fullFile - path to save the settings file
%
% RETURNS:
%   saves the settings of the simulation as a .txt file
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
    geometry    = model.experiment.geometry;
    rock        = model.experiment.rock;
    fluid       = model.experiment.fluid;
    process     = model.experiment.process;
    simulation  = model.simulation;
    schedule    = model.experiment.schedule;
    observation = model.experiment.observation;
    satfun      = model.experiment.satfun;
    plot        = model.plot;
    output      = model.output;
    history_match = model.history_match;
    
    %% file header
    file  = fopen(fullFile,'w');
    fprintf(file,'%s\n','# Comments starts with hash(#) symbol.');
    fprintf(file,'%s\n','# Between the keyword and the value should be at least one tab.');
    fprintf(file,'%s\n','# Press enter key at the end of each line of data, even the last line of the file.');
    fprintf(file,'%s\n','# Each block is separated by a blank line.');
    fprintf(file,'\n');
    
    %% block GEOMETRY = ["LENGTH" "DIAMETER"]
    blockHeader = 'GEOMETRY';
    WriteBlockHeader(file,blockHeader);
    if(isfield(geometry,'length'))
        if(geometry.length.include)
            fprintf(file,'%s','LENGTH');
            fprintf(file,'\t\t %g',geometry.length.inputValue);
            fprintf(file,'\t %s\n',geometry.length.inputUnit);
        end
    end
    if(isfield(geometry,'diameter'))
        if(geometry.diameter.include)
            fprintf(file,'%s','DIAMETER');
            fprintf(file,'\t %g',geometry.diameter.inputValue);
            fprintf(file,'\t %s\n',geometry.diameter.inputUnit);
        end
    end
    fprintf(file,'\n');
    
    %% block ROCK = ["POROSITY" "PERMEABILITY"]
    blockHeader = 'ROCK';
    WriteBlockHeader(file,blockHeader);
    if(isfield(rock,'poro'))
        if(rock.poro.include)
            fprintf(file,'%s','POROSITY');
            fprintf(file,'\t\t %g\n',rock.poro.inputValue);
        end
    end
    if(isfield(rock,'perm'))
        if(rock.perm.include)
            fprintf(file,'%s','PERMEABILITY');
            fprintf(file,'\t %g',rock.perm.inputValue);
            fprintf(file,'\t %s\n',rock.perm.inputUnit);
        end
    end
    if(isfield(rock,'Swi'))
        if(rock.Swi.include)
            fprintf(file,'%s','INITIALWATER');
            fprintf(file,'\t %g\n',rock.Swi.inputValue);
        end
    end
    fprintf(file,'\n');
    
    %% block FLUID = ["DENSITYW" "DENSITYNW" "VISCOSITYW" "VISCOSITYNW"]
    blockHeader = 'FLUID';
    WriteBlockHeader(file,blockHeader);
    if(isfield(fluid,'rhoW'))
        if(fluid.rhoW.include)
            fprintf(file,'%s','DENSITYW');
            fprintf(file,'\t %g',fluid.rhoW.inputValue);
            fprintf(file,'\t %s\n',fluid.rhoW.inputUnit);
        end
    end
    if(isfield(fluid,'rhoNW'))
        if(fluid.rhoNW.include)
            fprintf(file,'%s','DENSITYNW');
            fprintf(file,'\t %g',fluid.rhoNW.inputValue);
            fprintf(file,'\t %s\n',fluid.rhoNW.inputUnit);
        end
    end
    if(isfield(fluid,'muW'))
        if(fluid.muW.include)
            fprintf(file,'%s','VISCOSITYW');
            fprintf(file,'\t %g',fluid.muW.inputValue);
            fprintf(file,'\t\t %s\n',fluid.muW.inputUnit);
        end
    end
    if(isfield(fluid,'muNW'))
        if(fluid.muNW.include)
            fprintf(file,'%s','VISCOSITYNW');
            fprintf(file,'\t %g',fluid.muNW.inputValue);
            fprintf(file,'\t\t %s\n',fluid.muNW.inputUnit);
        end
    end
    fprintf(file,'\n');
    
    %% block PROCESS = ["SS" "USS" "CENT"]
    blockHeader = 'PROCESS';
    WriteBlockHeader(file,blockHeader);
    fprintf(file,'%s',process.type);
    fprintf(file,'\t\t %s\n',process.name);
    fprintf(file,'\n');
    
    %% block SIMULATION = ["TYPE" "NCELLS" "BCELLS" "TIMESTEP" "GAUGEOFF"]
    blockHeader = 'SIMULATION';
    WriteBlockHeader(file,blockHeader);
    fprintf(file,'%s','TYPE');
    fprintf(file,'\t\t %s\n',simulation.type);
    if(isfield(simulation,'nCells'))
        if(simulation.nCells.include)
            fprintf(file,'%s','NCELLS');
            fprintf(file,'\t\t %g\n',simulation.nCells.inputValue);
        end
    end
    if(isfield(simulation,'bCells'))
        if(simulation.bCells.include)
            fprintf(file,'%s','BCELLS');
            fprintf(file,'\t\t %g',simulation.bCells.inputValue);
            fprintf(file,'\t\t %g',simulation.bCells.firstCellSw.inputValue);
            fprintf(file,'\t\t %g\n',simulation.bCells.lastCellSw.inputValue);
        end
    end
    if(isfield(simulation,'timeStep'))
        if(simulation.timeStep.include)
            fprintf(file,'%s','MAXTIMESTEP');
            fprintf(file,'\t %g',simulation.timeStep.inputValue);
            fprintf(file,'\t %s\n',simulation.timeStep.inputUnit);
        end
    end
    if(isfield(simulation,'rampupsteps'))
        if(simulation.rampupsteps.include)
            fprintf(file,'%s','RAMPUPSTEPS');
            fprintf(file,'\t %g\n',simulation.rampupsteps.inputValue);
        end
    end
    if(isfield(simulation,'gaugeOff'))
        if(simulation.gaugeOff.include)
            fprintf(file,'%s','GAUGEOFF');
            fprintf(file,'\t %g',simulation.gaugeOff.inputValue);
            fprintf(file,'\t %s\n',simulation.gaugeOff.inputUnit);
        end
    end
    if(isfield(simulation,'high_percision_mode'))
        if simulation.high_percision_mode
            fprintf(file,'%s\n','HIGH_PERCISION_MODE');
        end
    end
            
    fprintf(file,'\n');

    %% block SCHEDULE = ["FILENAME" "PINIT" "POUT" "CENTRAD" "STARTUP"]
    blockHeader = 'SCHEDULE';
    WriteBlockHeader(file,blockHeader);
    if(isfield(schedule,'fileName'))
        fprintf(file,'%s','FILENAME');
        fprintf(file,'\t %s\n',fullfile(schedule.filePath,schedule.fileName));
    end
    if(isfield(schedule,'pini'))
        if(schedule.pini.include)
            fprintf(file,'%s','PINI');
            fprintf(file,'\t\t %g',schedule.pini.inputValue);
            fprintf(file,'\t\t %s\n',schedule.pini.inputUnit);
        end
    end
    if(isfield(schedule,'pout'))
        if(schedule.pout.include)
            fprintf(file,'%s','POUT');
            fprintf(file,'\t\t %g',schedule.pout.inputValue);
            fprintf(file,'\t\t %s\n',schedule.pout.inputUnit);
        end
    end
    if(strcmpi(process.type,'cent'))
        if isfield(schedule,'startup')
        if isfield(schedule.startup,'include') 
        if(schedule.startup.include)
            fprintf(file,'%s','STARTUP');
            fprintf(file,'\t\t %g',schedule.startupPeriod.inputValue);
            fprintf(file,'\t\t %s',schedule.startupPeriod.inputUnit);
            fprintf(file,'\t\t %g',schedule.startupOmega.inputValue);
            fprintf(file,'\t\t %s\n',schedule.startupOmega.inputUnit);
        end
        end
        end
        fprintf(file,'%s','CENTRAD');
        fprintf(file,'\t\t %g',schedule.centRad.inputValue);
        fprintf(file,'\t\t %s\n',schedule.centRad.inputUnit);
    end
    fprintf(file,'\n');    
    
    %% block OBSERVATION = ["PRESSURE" "SWAVG" "satProfile" "PRODUCTION"]
    blockHeader = 'OBSERVATION';
    WriteBlockHeader(file,blockHeader);
    if(isfield(observation,'pressure'))
        if(observation.pressure.include)
            fprintf(file,'%s','PRESSURE');
            fprintf(file,'\t\t %s\n',fullfile(observation.pressure.filePath,observation.pressure.fileName));
        end
    end
    if(isfield(observation,'swavg'))
        if(observation.swavg.include)
            fprintf(file,'%s','SWAVG');
            fprintf(file,'\t\t\t %s\n',fullfile(observation.swavg.filePath,observation.swavg.fileName));
        end
    end
    if(isfield(observation,'satProfile'))
        if(observation.satProfile.include)
            fprintf(file,'%s','SATPROFILE');
            fprintf(file,'\t\t %s\n',fullfile(observation.satProfile.filePath,observation.satProfile.fileName));
        end
    end
    if(isfield(observation,'prod'))
        if(observation.prod.include)
            fprintf(file,'%s','PRODUCTION');
            fprintf(file,'\t\t %s\n',fullfile(observation.prod.filePath,observation.prod.fileName));
        end
    end
    fprintf(file,'\n'); 
    
    %% block SATURATION FUNCTIONS = ["PC" "KRW" "KRNW" "SWC" "SOR" "KRW@SOR" "KRO@SWC"]
    blockHeader = 'SATURATION FUNCTIONS';
    WriteBlockHeader(file,blockHeader);
    kr = satfun.kr;
    if(kr.include)
        fprintf(file,'%s \t\t %s','KR',kr.type);
        if(strcmp(kr.type,strcat('TABLE')))
            fprintf(file,'\t\t %s\n',fullfile(kr.filePath,kr.fileName));
        end
        if(strcmp(kr.type,strcat('MODIFIED-COREY'))) 
            fprintf(file,'\t\t %g',kr.Swc);
            fprintf(file,'\t\t %g',kr.Sor);
            fprintf(file,'\t\t %g',kr.krwSor);
            fprintf(file,'\t\t %g',kr.kroSwc);     
            fprintf(file,'\t\t %g',kr.nW);
            fprintf(file,'\t\t %g\n',kr.nNW);     
        end
        if(strcmp(kr.type,strcat('BROOKS-COREY')))
            fprintf(file,'\t\t %g',kr.Swc);
            fprintf(file,'\t\t %g',kr.Sor);
            fprintf(file,'\t\t %g',kr.krwSor);
            fprintf(file,'\t\t %g',kr.kroSwc);
            fprintf(file,'\t\t %g\n',kr.lambda);
        end
        if(strcmp(kr.type,strcat('BURDINE')))
            fprintf(file,'\t\t %g',kr.krwSor);
            fprintf(file,'\t\t %g\n',kr.kroSwc);                   
        end
        if(strcmp(kr.type,strcat('LET')))
            fprintf(file,'\t\t %g',kr.Swc);
            fprintf(file,'\t\t %g',kr.Sor);
            fprintf(file,'\t\t %g',kr.krwSor);
            fprintf(file,'\t\t %g',kr.kroSwc);            
            fprintf(file,'\t\t %g',kr.Lw);
            fprintf(file,'\t\t %g',kr.Lnw);
            fprintf(file,'\t\t %g',kr.Ew);
            fprintf(file,'\t\t %g',kr.Enw);
            fprintf(file,'\t\t %g',kr.Tw);
            fprintf(file,'\t\t %g\n',kr.Tnw); 
        end
    end
    pc = satfun.pc;
    if(pc.include)
        if strcmpi(pc.type, 'zero')
           fprintf(file,'%s \t\t %s\n','PC',pc.type); 
        else
            fprintf(file,'%s \t\t %s','PC',pc.type);
        end
        if(strcmp(pc.type,strcat('TABLE')))
            fprintf(file,'\t\t\t\t %s\n',fullfile(pc.filePath,pc.fileName));
        end
        if(strcmp(pc.type,strcat('BROOKS-COREY')))
            fprintf(file,'\t\t %g',pc.pd);
            fprintf(file,'\t\t %g\n',pc.lambda);
        end
        if(strcmp(pc.type,strcat('SKJAEVELAND')))
            fprintf(file,'\t\t %g',pc.cwi);
            fprintf(file,'\t\t %g',pc.coi);
            fprintf(file,'\t\t %g',pc.awi);
            fprintf(file,'\t\t %g\n',pc.aoi);
        end
        if(strcmp(pc.type,strcat('MODIFIED-SKJAEVELAND')))
            fprintf(file,'\t\t %g',pc.cwi);
            fprintf(file,'\t\t %g',pc.coi);
            fprintf(file,'\t\t %g',pc.ri);
            fprintf(file,'\t\t %g',pc.bi);
            fprintf(file,'\t\t %g',pc.Swd);
            fprintf(file,'\t\t %g\n',pc.Sod);
        end
        if(strcmp(pc.type,strcat('LET-DRAINAGE')))
            fprintf(file,'\t\t %g',pc.max_pc);
            fprintf(file,'\t\t %g',pc.entry_pc);
            fprintf(file,'\t\t %g',pc.entry_multiplier);
            fprintf(file,'\t\t %g',pc.forced_multiplier);
            fprintf(file,'\t\t %g',pc.L_entry);
            fprintf(file,'\t\t %g',pc.E_entry);
            fprintf(file,'\t\t %g',pc.T_entry);
            fprintf(file,'\t\t %g',pc.L_forced);
            fprintf(file,'\t\t %g',pc.E_forced);
            fprintf(file,'\t\t %g\n',pc.T_forced);
        end
        if(strcmp(pc.type,strcat('LET-IMBIBITION')))
            fprintf(file,'\t\t %g',pc.sw_pc0);
            fprintf(file,'\t\t %g',pc.max_pc);
            fprintf(file,'\t\t %g',pc.min_pc);
            fprintf(file,'\t\t %g',pc.spontaneous_multiplier);
            fprintf(file,'\t\t %g',pc.forced_multiplier);
            fprintf(file,'\t\t %g',pc.L_spont);
            fprintf(file,'\t\t %g',pc.E_spont);
            fprintf(file,'\t\t %g',pc.T_spont);
            fprintf(file,'\t\t %g',pc.L_forced);
            fprintf(file,'\t\t %g',pc.E_forced);
            fprintf(file,'\t\t %g\n',pc.T_forced);
        end
    end
    fprintf(file,'\n'); 
    
    %% block PLOT OPTIONS = ["STYLE" "COLORMAP" "DISPLAYTIME"
    %                        "DISPLAYLENGTH" "DISPLAYVOLUME" 
    %                        "DISPLAYPRESS" "DISPLAYRATE"]
    blockHeader = 'PLOT OPTIONS';
    WriteBlockHeader(file,blockHeader);
    if(isfield(plot,'style'))
        if(plot.style.include)
            fprintf(file,'%s \t\t\t %s\n','STYLE',plot.style.inputStyle);
        end
    end
    if(isfield(plot,'colormap'))
        if(plot.colormap.include)
            fprintf(file,'%s \t\t %s\n','COLORMAP',plot.colormap.inputColormap);
        end
    end
    if(isfield(plot,'displaySat'))
        if(plot.displaySat.include)
            fprintf(file,'%s \t %s\n','DISPLAYSAT',plot.displaySat.inputUnit);
        end
    end
    if(isfield(plot,'displayTime'))
        if(plot.displayTime.include)
            fprintf(file,'%s \t %s\n','DISPLAYTIME',plot.displayTime.inputUnit);
        end
    end
    if(isfield(plot,'displayLength'))
        if(plot.displayLength.include)
            fprintf(file,'%s \t %s\n','DISPLAYLENGTH',plot.displayLength.inputUnit);
        end
    end
    if(isfield(plot,'displayVolume'))
        if(plot.displayVolume.include)
            fprintf(file,'%s \t %s\n','DISPLAYVOLUME',plot.displayVolume.inputUnit);
        end
    end
    if(isfield(plot,'displayPress'))
        if(plot.displayPress.include)
            fprintf(file,'%s \t %s\n','DISPLAYPRESS',plot.displayPress.inputUnit);
        end
    end
    if(isfield(plot,'displayRate'))
        if(plot.displayRate.include)
            fprintf(file,'%s \t %s\n','DISPLAYRATE',plot.displayRate.inputUnit);
        end
    end
    fprintf(file,'\n'); 
    
    %% block OUTPUT OPTIONS = ["SATPROFILE" "SAVECONFIG" "QUANTITIES"
    %                          "TIME" "SWAVG" "INJ" "PROD" "DELTAP"]
    blockHeader = 'OUTPUT OPTIONS';
    WriteBlockHeader(file,blockHeader);    
    if(isfield(output,'satProfile'))
        if(output.satProfile.include)
            fprintf(file,'%s \t\t %s\n','SATPROFILE',fullfile(output.satProfile.filePath,output.satProfile.fileName));
        end
    end   
    if(isfield(output,'saveConfig'))
        if(output.saveConfig.include)
            fprintf(file,'%s \t\t %s\n','SAVECONFIG',fullfile(output.saveConfig.filePath,output.saveConfig.fileName));
        end
    end
    if(isfield(output,'quantities')) && (output.quantities.include) && ...
       ( isfield(output.quantities,'filePath') && isfield(output.quantities,'fileName') )
        fprintf(file,'%s \t\t %s\n','QUANTITIES',fullfile(output.quantities.filePath,output.quantities.fileName));
    else
        fprintf(file,'%s \t\t %s\n','QUANTITIES','Unknown');    
    end
    if(isfield(output,'time'))
        if(output.time.include)
            fprintf(file,'%s\n','TIME');
        end
    end
    if(isfield(output,'swavg'))
        if(output.swavg.include)
            fprintf(file,'%s\n','SWAVG');
        end
    end
    if(isfield(output,'inj'))
        if(output.inj.include)
            fprintf(file,'%s\n','INJ');
        end
    end
    if(isfield(output,'prod'))
        if(output.prod.include)
            fprintf(file,'%s\n','PROD');
        end
    end
    if(isfield(output,'deltaP'))
        if(output.deltaP.include)
            fprintf(file,'%s\n','DELTAP');
        end
    end
    fprintf(file,'\n'); 
    
    %% block OBJECTIVE FUNCTION CONFIGURATIONS = ["EXCEL_FILE_NAME","EXCEL_FILE_PATH",
    %"MULTIPOINT","MULTIPOINT_RANDOM_NUMBERS","KR","KR_MODEL","PC","PC_MODEL",
    %"PDIFF_WEIGHT","SWAVG_WEIGHT","PROD_WEIGHT","SAT_PROFILE_WEIGHT","PDIF_ERROR",
    %"SWAVG_ERROR","PROD_ERROR","SAT_PROFILE_ERROR"]
    blockHeader = 'OBJECTIVE FUNCTION CONFIGURATIONS';
    WriteBlockHeader(file,blockHeader);    
    if(isfield(history_match,'hm_template_name'))
        fprintf(file,'%s \t\t %s\n','EXCEL_FILE_NAME',history_match.hm_template_name);
    end   
    if(isfield(history_match,'hm_template_path'))
        fprintf(file,'%s \t\t %s\n','EXCEL_FILE_PATH',history_match.hm_template_path);
    end   
    if(isfield(history_match,'multi_point'))
        if(history_match.multi_point)
            fprintf(file,'%s\n','MULTIPOINT');
        end
    end
    if(isfield(history_match,'kr'))
        if(history_match.kr.status)
            fprintf(file,'%s\n','KR');
        end
    end
    if(isfield(history_match,'kr'))
        if(history_match.kr.status)
            fprintf(file,'%s \t\t %s\n','KR_MODEL',history_match.kr.type);
        end
    end
    if(isfield(history_match,'pc'))
        if(history_match.pc.status)
            fprintf(file,'%s\n','PC');
        end
    end 
    if(isfield(history_match,'pc'))
        if(history_match.pc.status)
            fprintf(file,'%s \t\t %s\n','PC_MODEL',history_match.pc.type);
        end
    end    
    if(isfield(history_match,'pdiff_weight'))
        fprintf(file,'%s','PDIFF_WEIGHT');
        fprintf(file,'\t\t %g\n',history_match.pdiff_weight);
    end    
    if(isfield(history_match,'swavg_weight'))
        fprintf(file,'%s','SWAVG_WEIGHT');
        fprintf(file,'\t\t %g\n',history_match.swavg_weight);
    end    
    if(isfield(history_match,'prod_weight'))
        fprintf(file,'%s','PROD_WEIGHT');
        fprintf(file,'\t\t %g\n',history_match.prod_weight);
    end    
    if(isfield(history_match,'sat_profile_weight'))
        fprintf(file,'%s','SAT_PROFILE_WEIGHT');
        fprintf(file,'\t\t %g\n',history_match.sat_profile_weight);
    end    
    if(isfield(history_match,'pdiff_error'))
        fprintf(file,'%s','PDIF_ERROR');
        fprintf(file,'\t\t %g\n',history_match.pdiff_error);
    end    
    if(isfield(history_match,'swavg_error'))
        fprintf(file,'%s','SWAVG_ERROR');
        fprintf(file,'\t\t %g\n',history_match.swavg_error);
    end    
    if(isfield(history_match,'prod_error'))
        fprintf(file,'%s','PROD_ERROR');
        fprintf(file,'\t\t %g\n',history_match.prod_error);
    end    
    if(isfield(history_match,'sat_profile_error'))
        fprintf(file,'%s','SAT_PROFILE_ERROR');
        fprintf(file,'\t\t %g\n',history_match.sat_profile_error);
    end
    fprintf(file,'\n');
    
    %% block HISTORY MATCH CONFIGURATIONS = ["USE_PARALLEL","OPTIMALITY_TOLERANCE",
    %"STEP_TOLERANCE","MAX_FUNCTION_EVALUATIONS","SCALE_PROBLEM","OBJECTIVE_FUNCTION_TYPE",
    %"CENT_FILE_NAME","CENT_FILE_PATH","HISTORYMATCH_ALGORITHM"]
    
    blockHeader = 'HISTORY MATCH CONFIGURATIONS';
    WriteBlockHeader(file,blockHeader);    
    if(isfield(history_match,'UseParallel'))
        if history_match.UseParallel
            fprintf(file,'%s\n','USE_PARALLEL');
        end
    end       
    if(isfield(history_match,'OptimalityTolerance'))
        fprintf(file,'%s','OPTIMALITY_TOLERANCE');
        fprintf(file,'\t\t %g\n',history_match.OptimalityTolerance);
    end
    if(isfield(history_match,'StepTolerance'))
        fprintf(file,'%s','STEP_TOLERANCE');
        fprintf(file,'\t\t %g\n',history_match.StepTolerance);
    end        
    if(isfield(history_match,'MaxFunctionEvaluations'))
        fprintf(file,'%s','MAX_FUNCTION_EVALUATIONS');
        fprintf(file,'\t\t %g\n',history_match.MaxFunctionEvaluations);
    end        
    if(isfield(history_match,'ScaleProblem'))
        if history_match.ScaleProblem
            fprintf(file,'%s\n','SCALE_PROBLEM');
        end
    end        
    if(isfield(history_match,'obj_fun'))
        fprintf(file,'%s \t\t %s\n','OBJECTIVE_FUNCTION_TYPE',history_match.obj_fun);
    end        
    if(isfield(history_match,'Cent_file_name'))
        fprintf(file,'%s \t\t %s\n','CENT_FILE_NAME',history_match.Cent_file_name);
    end        
    if(isfield(history_match,'Cent_file_path'))
        fprintf(file,'%s \t\t %s\n','CENT_FILE_PATH',history_match.Cent_file_path);
    end    
    if(isfield(history_match,'algorithm'))
        fprintf(file,'%s \t\t %s\n','HISTORYMATCH_ALGORITHM',history_match.algorithm);
    end        
    fprintf(file,'\n');

    %% block MCMC SPECIFIC CONFIGURATIONS = ["RANDOM_SEED",
    %"SAMPLE_REFINEMENT_COUNT","CHAINSIZE","MPI_ENABLED"]    
    blockHeader = 'MCMC SPECIFIC CONFIGURATIONS';
    WriteBlockHeader(file,blockHeader);    
    if(isfield(history_match,'random_seed'))
        fprintf(file,'%s','RANDOM_SEED');
        fprintf(file,'\t\t %g\n',history_match.random_seed);
    end
    if(isfield(history_match,'sampleRefinementCount'))
        fprintf(file,'%s','SAMPLE_REFINEMENT_COUNT');
        fprintf(file,'\t\t %g\n',history_match.sampleRefinementCount);
    end
    if(isfield(history_match,'chainSize'))
        fprintf(file,'%s','CHAINSIZE');
        fprintf(file,'\t\t %g\n',history_match.chainSize);
    end
    if(isfield(history_match,'mpi_enabled'))
        if history_match.mpi_enabled
            fprintf(file,'%s\n','MPI_ENABLED');
        end
    end      
    fprintf(file,'\n');
    
    
    %% close the file
    fclose(file);        
end

function WriteBlockHeader(file,blockHeader)
    idy = 50;
    symbol = '-';
    line1 = repelem(symbol,idy);
    fprintf(file,'%s','#');
    fprintf(file,'%s\n',line1);
    fprintf(file,'%s','#');
    fprintf(file,'%s %s\n','',blockHeader);
    fprintf(file,'%s','#');
    fprintf(file,'%s\n',line1);
end