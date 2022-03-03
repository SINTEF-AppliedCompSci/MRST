%% Automated History Matching of SCAL Experiments
% 
%  Description : 
%
%  Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
%  History :
%  \change{1.0}{09-Nov-2021}{Original}
%
%  --------------------------------------------------
%  (c) 2021, Siroos Azizmohammadi,
%  Omidreza Amrollahinasab
%  Chair of Reservoir Engineering, University of Leoben, Austria
%  email: siroos.azizmohammadi@unileoben.ac.at
%  url: https://dpe.ac.at
%  --------------------------------------------------
%
%%
%% clear memory, close all figures and screen
clear; close all; clc

%% add MRST modules to working path
year = 2020; release = 'a';
mrstVersion = strcat(string(year),release);
main_script_dir = AddMRST(mrstVersion);
settings_dir = "W:\CO2_Displacement_paper_from_Holger\Decane_brine";

%% configure model from file
if exist('settings_dir','var')
    model = Configure(settings_dir,"settings_decane_brine2.txt");
else
    model = Configure(main_script_dir,"settings_EST_047_SS_Drain.txt");
end

%% App function
model.App.include = false;

%% verbose to display internal messages
verbose = false; 
model.verbose = verbose;

%% discretize geometry
model = CreateGrid(model);

%% create rock model
model = CreateRock(model);

%% create saturation functions model
model = CreatePc(model);
model = CreateKr(model);

%% plot satuartion functions
if not(strcmpi(model.simulation.type,strcat('historymatch')))
    PlotKrPc(model);
    if or(isfield(model.experiment.satfun,'kr_compare_1'), ...
            isfield(model.experiment.satfun,'pc_compare_1'))
        PlotKrPcComparison(model)
    end
    choice_1 = input('Check input saturation functions, continue?','s');
    if strcmpi(choice_1, 'no')
        return
    end
end

%% create fluid model
model = CreateFluid(model);

%% run Froward/Histormatch simulation
if not((strcmpi(model.simulation.type,strcat('historymatch'))))  
% ------------forward heterogenous model-----------------------------------
    if model.experiment.rock.heterogeneous
        fprintf('Simulating the prior homogeneous model...\n')
        model = Run(model, false);
        fprintf('calculating the pc scaling factors...\n')
        alpha = model.experiment.rock.alpha.value;
        t_obs = model.experiment.observation.satProfile.table(2:end,1);
        index_mask = model.experiment.rock.het_index_mask;
        index_mask_time = t_obs(index_mask) / 3600;
        fmt = ['You are using the saturation profiles at ' repmat(' %.2f ',...
            1,numel(index_mask_time)) ' (hours) for f factor calculations\n'];
        fprintf(fmt, index_mask_time)
        f_factor = f_factor_calculator(model, alpha, true, model.experiment.rock.het_index_mask);
        choice_2 = input("Check f factor calculations, continue?", "s");
        if strcmpi(choice_2, "no")
            return
        end
        model = CreateGrid_heterogeneous(model);
        model = CreateRock_heterogeneous(model, f_factor);
        model = CreateFluid_heterogenous(model, f_factor);
        fprintf('Simulating the heterogeneous model...\n')
        model = Run(model, true);
        choice_3 = input("Plot saturation profile?", "s");
        if strcmpi(choice_3, "yes")
            plot_saturation_profile_modified(model, model.experiment.rock.het_index_mask)
        end
    else
% ------------forward homogeneous model-----------------------------------
        fprintf('Simulating the homogeneous model...\n')
        model = Run(model, true);
        choice_4 = input("Plot saturation profile?", "s");
        if strcmpi(choice_4, "yes")
            plot_saturation_profile(model)
        end
    end
else
% ------------investigate the input range for------------------------------
%     plot_input_range(model)
%     fprintf('Want to continue? Press Enter to continue \n')
%     pause
%-------------Exploratory data analysis------------------------------------
    filter_no.pressure = 50;
    filter_no.prod = 10;
    filtered_data = PlotObservation_pre_hm(model, filter_no);
    plot_saturation_profile_pre_hm(model)
    choice = input('History match with the filtered data?','s');
    if strcmpi(choice, 'yes')
        model = replace_observation_with_filtered_data(model, filtered_data);
    end
    close
%--------------------------------------------------------------------------
    fprintf('Starting the history match...\n')
%----------Multi objective history matching ------------------------------
    if strcmpi(model.history_match.algorithm,'ga_multi_objective')
        model = historymatch_multi_obj(model);
        [x_best, fval_sum_best] = get_x_best(model);
        model.history_match.x = x_best;
        model.history_match.fval = fval_sum_best;
    else
%----------single objective history matching ------------------------------
    [model, model_cent] = historymatch(model);
    end
    plot_HM_results(model)
    obj_fun = model.history_match.obj_fun;
    if strcmp(obj_fun,'Simultaneous')
        plot_HM_response(model,'ss')
        model_cent.history_match = model.history_match;
        plot_HM_response(model_cent,'cent')
    else
        if model.experiment.rock.heterogeneous
            plot_HM_response(model,'', model.experiment.rock.het_index_mask)
        else
            plot_HM_response(model,'')
        end
    end
end

%% output results
% if not(isempty(model.output))
%     if(model.output.include)
%         SaveResults(model);
%     end
% end

%% compute and plot fractional flow
% [sw, fw] = compute_fractional_flow (model);
% plot_fractional_flow(sw, fw);

