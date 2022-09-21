%% Automated History Matching of SCAL Experiments
%
% DESCRIPTION: Example of how to run a case with the SCAL module
% Have a look at "docs/list of examples.docx" for a list of available input
% examples + "doc/List of keywords in the settings files.docx" for a
% documentation about what is the role of each keyword
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

%% configure model from file
% adding mrst modules
mrstModule add ad-core ad-props ad-blackoil ad-scal
% input absolute path to the settings file containing the keywords and
% input parameters/tables
model = Configure(fullfile (ROOTDIR,'modules','ad-scal','examples','settings_case1.txt'));

%% verbose to display internal messages
% set true to show detailed messages about the simulation
verbose = false; 
model.verbose = verbose;

%% discretize geometry
% create the grid on which the simulation is run
model = CreateGrid(model);

%% create rock model
% setting porosity and permeability data for the simulation
model = CreateRock(model);

%% create saturation functions model
% creating the capillary pressure and relative permeability tables for the
% simulation
model = CreatePc(model);
model = CreateKr(model);

%% plot satuartion functions and comparison plots
% plot the pc and kr tables created in the previous step for verification
% and comparison purposes 
if not(strcmpi(model.simulation.type,strcat('historymatch')))
    PlotKrPc(model);
    if or(isfield(model.experiment.satfun,'kr_compare'), ...
            isfield(model.experiment.satfun,'pc_compare'))
        PlotKrPcComparison(model)
    end
    % just a pause before you run the simulations to check if you have the
    % correct input relative permeability and capillary pressures for 
    % your simulations - you can also plot them against other kr and pc 
    % curves for comparison purposes
    disp('Check input saturation functions, continue? [enter for yes/control + c to no]');
    pause()
end

%% create fluid model
% create functions to interpolate the pc and kr saturation function tables
model = CreateFluid(model);

%% run Froward/Histormatch simulation
% run the forward simulation: meaning that we have the pc and kr curve as
% input and we want to predict: pressure, production and saturation
% profiles
if not((strcmpi(model.simulation.type,strcat('historymatch'))))  
% ------------forward homogeneous model-----------------------------------
    fprintf('Simulating the homogeneous model...\n')
    model = Run(model, true);
    disp("Plot saturation profile? [enter for yes/ control + c for no]");
    pause()
    plot_saturation_profile(model)
else
% ------------investigate the input range for------------------------------
% for the history matching part, it is usefull to plot the range of kr and
% pc tables that we are considering as boundary conditions before running
% the history matching simulations
%     plot_input_range(model)
%     fprintf('Want to continue? Press Enter to continue \n')
%     pause
%-------------Exploratory data analysis------------------------------------
% in case of noisy data, there is the option to filter the noise out using
% a median filter, the higher the number, the more aggressive the filtering
% will be, and then the history match can be done with the filtered data or
% the original one
    filter_no.pressure = 5;
    filter_no.prod = 5;
    filtered_data = PlotObservation_pre_hm(model, filter_no);
    plot_saturation_profile_pre_hm(model)
    choice = input('History match with the filtered data? [y/n]','s');
    if strcmpi(choice, 'y')
        model = replace_observation_with_filtered_data(model, filtered_data);
    end
    close
%--------------------------------------------------------------------------
    fprintf('Starting the history match...\n')
% for the history matching part, it can be done using multi objective
% genetic algorithm or single objective algorithms like active-set in which
% we sum up the errors from different measurements/parameters
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
% visualize the results of history match simulations and compare with 
% measurements  
    plot_HM_results(model)
    obj_fun = model.history_match.obj_fun;
    if strcmpi(obj_fun,'simultaneous')
        plot_HM_response(model,'ss')
        model_cent.history_match = model.history_match;
        plot_HM_response(model_cent,'cent')
    else
        plot_HM_response(model,'')
    end
end
%% compute and plot fractional flow
% plot and compute the fractional flow
% most useful for steady state experimetns 
% [sw, fw] = compute_fractional_flow (model);
% plot_fractional_flow(sw, fw);

