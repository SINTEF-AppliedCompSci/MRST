function PlotKrPcComparison(model)
%
% DESCRIPTION: compares the relative permeability used for the modeling and
%              the ones input with KR_COMPARE or PC_COMPARE keywords
%
% SYNOPSIS:
%   PlotKrPcComparison(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling and the 
%                 comparison tables
%
% RETURNS:
%   plot of the capillary pressure and relative permeability used for the
%   model and the comparisons
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
    satfun_original_model = model.satfun;

    satfun = sat_fun_creation(model);

    figure('Name', 'Saturation Functions comparison', ...
            'NumberTitle', 'off', ...
            'WindowStyle', 'docked' );

    tiledlayout("flow")
    font_size = 10;

    %% relative permeabilities (linear)
    nexttile; hold on
    kr_legend_array = {};
    p_kr_1 = plot(satfun_original_model.sw_kr, satfun_original_model.krw ,'r-');
    kr_legend_array{end + 1} = 'kr model';
    p_kr_2 = plot(satfun_original_model.sw_kr, satfun_original_model.kro ,'r-');
    kr_legend_array{end + 1} = '';
    % we use this one to fix double legends
    kr_legend_counter = 0;
    % to fill the space between upper and lower boundary
    kr_coloring_counter = 0;
    % to use another plotting config when several sat functions are compared
    compare_counter = 0;
    for i = 1:length(satfun)
        if strcmpi(satfun{i}.name, 'kr_compare')
            if compare_counter == 0
                p_kr_3 = plot(satfun{i}.sw_kr, satfun{i}.krw ,'ks');
                kr_legend_array{end + 1} = 'kr compare';
                p_kr_4 = plot(satfun{i}.sw_kr, satfun{i}.kro ,'ks');
                kr_legend_array{end + 1} = '';
            else
                p_kr_3 = plot(satfun{i}.sw_kr, satfun{i}.krw ,'b-');
                kr_legend_array{end + 1} = 'kr compare';
                p_kr_4 = plot(satfun{i}.sw_kr, satfun{i}.kro ,'b-');
                kr_legend_array{end + 1} = '';
            end
            compare_counter = compare_counter + 1;
        end
        if strcmpi(satfun{i}.name, 'kr_boundary')
            p_kr_5 = plot(satfun{i}.sw_kr, satfun{i}.krw ,'k--');
            if kr_legend_counter == 0
                kr_legend_array{end + 1} = 'Uncertainty boundary';
                kr_legend_counter = 1;
            else
                kr_legend_array{end + 1} = '';
            end
            p_kr_6 = plot(satfun{i}.sw_kr, satfun{i}.kro ,'k--');
            kr_legend_array{end + 1} = '';
            kr_coloring_counter = kr_coloring_counter + 1;
            if rem(kr_coloring_counter,2) == 0
                xx = [satfun{i}.sw_kr', fliplr(satfun{i}.sw_kr')];
                yy_krw = [ykrw_old', fliplr(satfun{i}.krw')];
                yy_kro = [ykro_old', fliplr(satfun{i}.kro')];
                fill(xx, yy_krw, 'b','FaceAlpha',0.1)
                fill(xx, yy_kro, 'g','FaceAlpha',0.1)
            else
                ykrw_old = satfun{i}.krw;
                ykro_old = satfun{i}.kro;
            end
        end
    end
    legend(kr_legend_array{:})
    legend('Location','best','FontSize',font_size); legend boxoff
    set(gca,'FontSize', font_size)
    xlabel('Water saturation')
    ylabel('Relative permeability')
    xlim([0 1])
    ylim([0 1])
    grid off

    %% relative permeabilities (logarithmic)
    nexttile; hold on
    kr_legend_array = {};
    p_kr_1 = plot(satfun_original_model.sw_kr, satfun_original_model.krw ,'r-');
    kr_legend_array{end + 1} = 'kr model';
    p_kr_2 = plot(satfun_original_model.sw_kr, satfun_original_model.kro ,'r-');
    kr_legend_array{end + 1} = '';
    kr_legend_counter = 0;
    kr_coloring_counter = 0;
    compare_counter = 0;
    for i = 1:length(satfun)
        mask = satfun{i}.krw > 0 & satfun{i}.kro > 0;
        sw_kr = satfun{i}.sw_kr(mask);
        krw = satfun{i}.krw(mask);
        kro = satfun{i}.kro(mask);
        if strcmpi(satfun{i}.name, 'kr_compare')
            if compare_counter == 0
                p_kr_3 = plot(sw_kr, krw ,'ks');
                kr_legend_array{end + 1} = 'kr compare';
                p_kr_4 = plot(sw_kr, kro ,'ks');
                kr_legend_array{end + 1} = '';
            else
                p_kr_3 = plot(sw_kr, krw ,'b-');
                kr_legend_array{end + 1} = 'kr compare';
                p_kr_4 = plot(sw_kr, kro ,'b-');
                kr_legend_array{end + 1} = '';
            end
            compare_counter = compare_counter + 1;
        end
        if strcmpi(satfun{i}.name, 'kr_boundary')
            p_kr_5 = plot(sw_kr, krw ,'k--');
            if kr_legend_counter == 0
                kr_legend_array{end + 1} = 'Uncertainty boundary';
                kr_legend_counter = 1;
            else
                kr_legend_array{end + 1} = '';
            end
            p_kr_6 = plot(sw_kr, kro ,'k--');
            kr_legend_array{end + 1} = '';
            kr_coloring_counter = kr_coloring_counter + 1;
            if rem(kr_coloring_counter,2) == 0
                xx = [sw_kr_old', fliplr(sw_kr')];
                yy_krw = [ykrw_old', fliplr(krw')];
                yy_kro = [ykro_old', fliplr(kro')];
                fill(xx, yy_krw, 'b','FaceAlpha',0.1)
                fill(xx, yy_kro, 'g','FaceAlpha',0.1)
            else
                sw_kr_old = sw_kr;
                ykrw_old = krw;
                ykro_old = kro;
            end
        end
    end
    legend(kr_legend_array{:})
    legend('Location','best','FontSize',font_size); legend boxoff
    set(gca,'FontSize', font_size)
    xlabel('Water saturation')
    ylabel('Relative permeability')
    xlim([0 1])
    ylim([1e-5 1])
    grid off
    set(gca, 'YScale', 'log')

    %% capillary pressure (linear)
    nexttile; hold on
    pc_legend_array = {};
    p_pc_1 = plot(satfun_original_model.sw_pc, satfun_original_model.pc/...
        Convert('bar'),'r-');
    pc_legend_array{end + 1} = 'pc model';
    pc_boundary_counter = 0;
    pc_coloring_counter = 0;
    for i = 1:length(satfun)
        if strcmpi(satfun{i}.name, 'pc_compare')
            p_pc_2 = plot(satfun{i}.sw_pc, satfun{i}.pc /...
                Convert('bar'),'ks');
            pc_legend_array{end + 1} = 'pc compare';
        end
        if strcmpi(satfun{i}.name, 'pc_boundary')
            p_pc_3 = plot(satfun{i}.sw_pc, satfun{i}.pc / ...
                Convert('bar'),'k--');
            if pc_boundary_counter == 0
                pc_legend_array{end + 1} = 'Uncertainty boundary';
                pc_boundary_counter = 1;
            else
                pc_legend_array{end + 1} = '';
            end
            pc_coloring_counter = pc_coloring_counter + 1;
            if rem(pc_coloring_counter,2) == 0
                xx = [sw_pc_old', fliplr(satfun{i}.sw_pc')];
                yy_pc = [ypc_old', fliplr(satfun{i}.pc')];
                fill(xx, yy_pc / Convert('bar'), 'k','FaceAlpha',0.1)
            else
                sw_pc_old = satfun{i}.sw_pc;
                ypc_old = satfun{i}.pc;
            end
        end
    end
    legend(pc_legend_array{:})
    legend('Location','best','FontSize',font_size); legend boxoff
    xlabel('Water saturation')
    ylabel('Capillary pressure (bar)')
    xlim([0 1])
    data = satfun_original_model.pc(:) / Convert('bar');
    for i = 1:length(satfun)
        data = [data ;satfun{i}.pc(:)] / Convert('bar');
    end
    if min(data) >= 0
        ylim([min(data)-0.1 min(max(data),2) + 0.1]);
    else
        ylim([max(min(data),-2) min(max(data),2)]);
    end
    set(gca,'FontSize', font_size)
    grid off

    %% capillary pressure (log)
    nexttile; hold on
    pc_legend_array = {};
    p_pc_1 = plot(satfun_original_model.sw_pc, satfun_original_model.pc/...
        Convert('bar'),'r-');
    pc_legend_array{end + 1} = 'pc model';
    pc_boundary_counter = 0;
    for i = 1:length(satfun)
        mask = satfun{i}.pc > 0;
        sw_pc = satfun{i}.sw_pc(mask);
        pc = satfun{i}.pc(mask);
        if strcmpi(satfun{i}.name, 'pc_compare')
            p_pc_2 = plot(sw_pc, pc /...
                Convert('bar'),'ks');
            pc_legend_array{end + 1} = 'pc compare';
        end
        if strcmpi(satfun{i}.name, 'pc_boundary')
            p_pc_3 = plot(sw_pc, pc / ...
                Convert('bar'),'k--');
            if pc_boundary_counter == 0
                pc_legend_array{end + 1} = 'Uncertainty boundary';
                pc_boundary_counter = 1;
            else
                pc_legend_array{end + 1} = '';
            end
            pc_coloring_counter = pc_coloring_counter + 1;
            if rem(pc_coloring_counter,2) == 0
                xx = [sw_pc_old', fliplr(sw_pc')];
                yy_pc = [ypc_old', fliplr(pc')];
                fill(xx, yy_pc / Convert('bar'), 'k','FaceAlpha',0.1)
            else
                sw_pc_old = sw_pc;
                ypc_old = pc;
            end
        end
    end
    legend(pc_legend_array{:})
    legend('Location','best','FontSize',font_size); legend boxoff
    xlabel('Water saturation')
    ylabel('Capillary pressure (bar)')
    xlim([0 1])
    data = satfun_original_model.pc(:) / Convert('bar');
    for i = 1:length(satfun)
        data = [data ;satfun{i}.pc(:)] / Convert('bar');
    end
    if min(data) >= 0
        ylim([0 min(max(data),2) + 0.1]);
    else
        ylim([0 min(max(data),2)]);
    end
    set(gca,'FontSize', font_size)
    grid off
    set(gca, 'YScale', 'log')


end

function satfun = sat_fun_creation(model)
%
% DESCRIPTION: creates the tables for the comparisons of relative
%              permeability and capillary pressure
%
% SYNOPSIS:
%   satfun = sat_fun_creation(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling and the 
%                 comparison tables
%
% RETURNS:
%   satfun - struct with relative permeability and capillary pressure
%   tables
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
    model_copy = model;
    satfun = {};
    if isfield(model.experiment.satfun, 'pc_compare')
        for i = 1:length(model.experiment.satfun.pc_compare)
            model_copy = model;
            model_copy.experiment.satfun.pc = model.experiment.satfun.pc_compare{i};
            model_copy = CreatePc(model_copy);
            satfun{end + 1} = model_copy.satfun;
            satfun{end}.name = 'pc_compare';
        end
    end
    if isfield(model.experiment.satfun, 'kr_compare')
        model_copy.experiment.satfun.kr = model.experiment.satfun.kr_compare;
        for i = 1:length(model.experiment.satfun.kr_compare)
            model_copy = model;
            model_copy.experiment.satfun.kr = model.experiment.satfun.kr_compare{i};
            model_copy = CreateKr(model_copy);
            satfun{end + 1} = model_copy.satfun;
            satfun{end}.name = 'kr_compare';
        end
    end
    if isfield(model.experiment.satfun, 'pc_boundary')
        model_copy.experiment.satfun.pc = model.experiment.satfun.pc_boundary;
        for i = 1:length(model.experiment.satfun.pc_boundary)
            model_copy = model;
            model_copy.experiment.satfun.pc = model.experiment.satfun.pc_boundary{i};
            model_copy = CreatePc(model_copy);
            satfun{end + 1} = model_copy.satfun;
            satfun{end}.name = 'pc_boundary';
        end
    end
    if isfield(model.experiment.satfun, 'kr_boundary')
        model_copy.experiment.satfun.kr = model.experiment.satfun.kr_boundary;
        for i = 1:length(model.experiment.satfun.kr_boundary)
            model_copy = model;
            model_copy.experiment.satfun.kr = model.experiment.satfun.kr_boundary{i};
            model_copy = CreateKr(model_copy);
            satfun{end + 1} = model_copy.satfun;
            satfun{end}.name = 'kr_boundary';
        end
    end
end
