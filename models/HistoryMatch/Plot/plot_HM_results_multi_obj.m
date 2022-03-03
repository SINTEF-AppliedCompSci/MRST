function plot_HM_results_multi_obj(model,index)
 
    sim_result = model.history_match.x(index,:); 
    rng default
    color = rand(height(sim_result),3);
    for i = 1:height(sim_result)
        x = sim_result(i,:);
        [sw_sim_result_kr,krw_sim_result,kro_sim_result,sw_sim_result_pc,sim_result_pc]  = saturation_tables(x,model);
        if i == 1
            figure('Name','Capillary Pressure','NumberTitle','off')
            ax_pc = gca;
            title('Capillary Pressure')
            hold on
            axis tight
            xlim([0 1])
            ylabel('Capillary Pressure (bar)')
            xlabel('Water Saturation (fraction)')
        end
        plot(ax_pc, sw_sim_result_pc,sim_result_pc*1e-5,'b-o', 'MarkerSize',2, 'color', color(i,:))

        if i == 1
            figure('Name','Relative Permeability','NumberTitle','off')
            ax_kr = gca;
            title('Relative Permeability')
            hold on
            ylim([0 1])
            xlim([0 1])
            ylabel('Relative Permeability')
            xlabel('Water Saturation (fraction)')
        end
        plot(ax_kr, sw_sim_result_kr,krw_sim_result, 'b-o', sw_sim_result_kr,kro_sim_result, 'b-o', 'MarkerSize',2, 'color', color(i,:))

    end
end

function [Sw,krw,kro,Sw_pc,pc] = saturation_tables(x,model)

    model = CreateSatFun_HM(x,model);
    Sw = model.satfun.sw_kr;
    Sw_pc = model.satfun.sw_pc;
    pc = model.satfun.pc;
    krw = model.satfun.krw;
    kro = model.satfun.kro;
end
