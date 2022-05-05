%% Automated History Matching of SCAL Experiments
%
% DESCRIPTION: Example of how to run a case with the SCAL module
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
%% clear memory, close all figures and screen
clc;
clear all;
close all;
clear classes;

%% add MRST modules to working path
year = 2020; release = 'a';
mrstVersion = strcat(string(year),release);
AddMRST(mrstVersion);
% in case it is in the main dir, user can use ".\"
% settings_dir = "W:\CO2_Displacement_paper_from_Holger\CO2_brine";
settings_dir = "W:\CO2_Displacement_paper_from_Holger\Decane_brine";
% settings_dir = "./";

%% configure model from file
model = Configure(settings_dir,"settings_decane_brine2.txt");

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

%% plot satuartion functions and comparison plots
if not(strcmpi(model.simulation.type,strcat('historymatch')))
    PlotKrPc(model);
    if or(isfield(model.experiment.satfun,'kr_compare'), ...
            isfield(model.experiment.satfun,'pc_compare'))
        PlotKrPcComparison(model)
    end
    choice_1 = input('Check input saturation functions, continue? [y/n]','s');
    if strcmpi(choice_1, 'n')
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
        choice_2 = input("Check f factor calculations, continue? [y/n]", "s");
        if strcmpi(choice_2, "n")
            return
        end
        model = CreateGrid_heterogeneous(model);
        model = CreateRock_heterogeneous(model, f_factor);
        model = CreateFluid_heterogenous(model, f_factor);
        fprintf('Simulating the heterogeneous model...\n')
        model = Run(model, true);
        choice_3 = input("Plot saturation profile? [y/n]", "s");
        if strcmpi(choice_3, "y")
            plot_saturation_profile_modified(model, model.experiment.rock.het_index_mask)
        end
    else
% ------------forward homogeneous model-----------------------------------
        fprintf('Simulating the homogeneous model...\n')
        model = Run(model, true);
        choice_4 = input("Plot saturation profile? [y/n]", "s");
        if strcmpi(choice_4, "y")
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
    choice = input('History match with the filtered data? [y/n]','s');
    if strcmpi(choice, 'y')
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
%% plotting after history match
    plot_HM_results(model)
    obj_fun = model.history_match.obj_fun;
    if strcmpi(obj_fun,'simultaneous')
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

