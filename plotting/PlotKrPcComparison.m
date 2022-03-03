function PlotKrPcComparison(model)

    satfun_original_model = model.satfun;

    satfun_compare_1 = sat_fun_creation(model, 'main');
    satfun_ub = sat_fun_creation(model, 'ub');
    satfun_lb = sat_fun_creation(model, 'lb');

    figure('Name', 'Saturation Functions comparison', ...
            'NumberTitle', 'off', ...
            'WindowStyle', 'docked' );

    tiledlayout(1,2)
    font_size = 16;
    nexttile; hold on
    kr_legend_array = [];
    p_kr_1 = plot(satfun_original_model.sw_kr, satfun_original_model.krw ,'r-');
    p_kr_2 = plot(satfun_original_model.sw_kr, satfun_original_model.kro ,'r-');
    kr_legend_array(end + 1) = p_kr_1;
    if isfield(model.experiment.satfun, 'kr_compare_1')
        p_kr_3 = plot(satfun_compare_1.sw_kr, satfun_compare_1.krw ,'ks');
        p_kr_4 = plot(satfun_compare_1.sw_kr, satfun_compare_1.kro ,'ks');
        kr_legend_array(end + 1) = p_kr_3;
    end
    if isfield(model.experiment.satfun, 'kr_ub')
        p_kr_5 = plot(satfun_ub.sw_kr, satfun_ub.krw ,'k--'); 
        p_kr_6 = plot(satfun_ub.sw_kr, satfun_ub.kro ,'k--');
        kr_legend_array(end + 1) = p_kr_5;
    end
    if isfield(model.experiment.satfun, 'kr_lb')
        p_kr_7 = plot(satfun_lb.sw_kr, satfun_lb.krw,'k--');
        p_kr_8 = plot(satfun_lb.sw_kr, satfun_lb.kro,'k--');
    end
    kr_label_cell = {'kr model', 'kr compare 1', 'Uncertainty boundary'};
    legend(kr_legend_array,kr_label_cell(1:length(kr_legend_array)));
    legend show; legend('Location','best','FontSize',font_size); legend boxoff
    set(gca,'FontSize', font_size)
    xlabel('Water saturation')
    ylabel('Relative permeability')
    xlim([0 1])
    ylim([0 1])
    grid off

    nexttile; hold on
    pc_legend_array = [];
    p_pc_1 = plot(satfun_original_model.sw_pc, satfun_original_model.pc / Convert('bar'),'r-');
    pc_legend_array(end + 1) = p_pc_1;
    if isfield(model.experiment.satfun, 'pc_compare_1')
        p_pc_2 = plot(satfun_compare_1.sw_pc, satfun_compare_1.pc / Convert('bar'),'ks');
        pc_legend_array(end + 1) = p_pc_2;
    end
    if isfield(model.experiment.satfun, 'pc_ub')
        p_pc_3 = plot(satfun_ub.sw_pc, satfun_ub.pc / Convert('bar'),'k--');
        pc_legend_array(end + 1) = p_pc_3;
    end
    if isfield(model.experiment.satfun, 'pc_lb')
        p_pc_4 = plot(satfun_lb.sw_pc, satfun_lb.pc / Convert('bar'),'k--');
    end
    pc_label_cell = {'pc model', 'pc compare 1', 'Uncertainty boundary'};
    legend(pc_legend_array, pc_label_cell(1:length(pc_legend_array)));
    legend show; legend('Location','best','FontSize',font_size); legend boxoff
    xlabel('Water saturation')
    ylabel('Capillary pressure (bar)')
    xlim([0 1])
    data = [satfun_original_model.pc(:);satfun_compare_1.pc(:)] / Convert('bar');
    if min(data) >= 0
        ylim([min(data)-0.1 min(max(data),2) + 0.1]);
    else
        ylim([max(min(data),-2) min(max(data),2)]);
    end
    set(gca,'FontSize', font_size)
    grid off
end

function satfun = sat_fun_creation(model, type)

    model_copy = model;
    if strcmpi(type, 'main')
        if isfield(model.experiment.satfun, 'pc_compare_1')
            model_copy.experiment.satfun.pc = model.experiment.satfun.pc_compare_1;
        end
        if isfield(model.experiment.satfun, 'kr_compare_1')
            model_copy.experiment.satfun.kr = model.experiment.satfun.kr_compare_1;
        end
    elseif strcmpi(type, 'ub')
        if isfield(model.experiment.satfun, 'pc_ub')
            model_copy.experiment.satfun.pc = model.experiment.satfun.pc_ub;
        end
        if isfield(model.experiment.satfun, 'kr_ub')
            model_copy.experiment.satfun.kr = model.experiment.satfun.kr_ub;
        end
    elseif strcmpi(type, 'lb')
        if isfield(model.experiment.satfun, 'pc_lb')
            model_copy.experiment.satfun.pc = model.experiment.satfun.pc_lb;
        end
        if isfield(model.experiment.satfun, 'kr_lb')
            model_copy.experiment.satfun.kr = model.experiment.satfun.kr_lb;
        end
    end
    model_copy = CreatePc(model_copy);
    model_copy = CreateKr(model_copy);
    satfun = model_copy.satfun;
end
