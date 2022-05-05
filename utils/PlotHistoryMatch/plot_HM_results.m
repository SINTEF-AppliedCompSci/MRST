function plot_HM_results(model)
%
% DESCRIPTION: after the history matching is done, plots the resulting
%              optimum saturation profiles
%
% SYNOPSIS:
%   plot_HM_results(model)
%
% PARAMETERS:
%   - model - struct after the history matching is done
%
% RETURNS:
%   plots the optimum saturation profile calculated with the history
%   matching algorithm
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
    sim_result = model.history_match.x; 
    x0 = model.history_match.x0;
     
    [sw_x0_kr,krw_x0,kro_x0,sw_x0_pc,x0_pc]  = saturation_tables(x0,model);
    [sw_hm_result_kr,krw_hm_result,kro_hm_result,sw_hm_result_pc,hm_result_pc]...
        = saturation_tables(sim_result,model);
    
    f = figure('Name','History matching results','NumberTitle','off');
    f.Position  = [500 500 1250 500];
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile;
    title('Capillary Pressure')
    hold on
    plot(sw_x0_pc, x0_pc,'r-o', 'MarkerSize',2)
    plot(sw_hm_result_pc, hm_result_pc,'b-o', 'MarkerSize',2)
    axis tight; grid on;
    xlim([0 1])
    if min(hm_result_pc) >= 0
        ylim([min(hm_result_pc)-0.1 min(max(hm_result_pc),2) + 0.1]);
    else
        ylim([min(min(hm_result_pc),-2) min(max(hm_result_pc),2)]);
    end
    legend('Initial Point', 'History match Result')
    ylabel('Capillary Pressure (bar)')
    xlabel('Water Saturation (fraction)')
 
    nexttile;
    title('Relative Permeability')
    hold on; grid on;
    plot(sw_x0_kr,krw_x0, 'r-o', sw_x0_kr,kro_x0, 'r-o', 'MarkerSize',2)
    plot(sw_hm_result_kr,krw_hm_result, 'b-o', sw_hm_result_kr,kro_hm_result, 'b-o', 'MarkerSize',2)
    ylim([0 1])
    xlim([0 1])
    legend('Initial Point', '', 'History match Result','')
    ylabel('Relative Permeability')
    xlabel('Water Saturation (fraction)')
end

function [Sw,krw,kro,Sw_pc,pc] = saturation_tables(x,model)

    params_no_kr = get_params_no_kr(model);
    model = Create_pc_history_match(x, model, params_no_kr);
    model = Create_kr_history_match(x, model);
    Sw = model.satfun.sw_kr;
    Sw_pc = model.satfun.sw_pc;
    pc = model.satfun.pc / Convert('bar');
    krw = model.satfun.krw;
    kro = model.satfun.kro;
end
