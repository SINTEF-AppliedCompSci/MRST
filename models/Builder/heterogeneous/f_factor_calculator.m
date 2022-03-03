function f_factor = f_factor_calculator(model, alpha, plots, index_mask)
% <keywords>
%
% Purpose : calculate the capillary pressure scaling factors
%
% Syntax : 
%   f_factor = f_factor_calculator(model)
%
% Input Parameters :
%   model: struct containing modeling params
%   plots: boolean to activate diagnostic plots
%
% Return Parameters :
%   f_factor: array containing capillary pressure scaling factors
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%


    % get the saturation from the observation
    x_obs = model.experiment.observation.satProfile.table(1,2:end);
    sw_obs = model.experiment.observation.satProfile.table(2:end,2:end);
    t_obs = model.experiment.observation.satProfile.table(2:end,1);
    
    % get the pressure from the simulation
    x = model.history_match.Sw_profile(1,2:end);
    t = model.history_match.Sw_profile(2:end,1);
    sw_model = model.history_match.Sw_profile(2:end,2:end);
    
    sw_obs_interped = interp2(x_obs, t_obs, sw_obs, x, t_obs, 'linear');
    sw_model_interped = interp2(x, t, sw_model, x, t_obs, 'linear');
    
% pressure from modeling can either come from pressure profile or using the
% capillary pressure curve (there is a slight difference between, reason 
% is unknown)
%     p = model.history_match.pressure_prof(2:end,2:end);
%     pout = model.experiment.schedule.pout.value;
%     dp_model = interp2(x, t, p, x_obs, t_obs, 'linear') - pout;
%     pc_model = dp_model;

    pc_model = model.fluid.pcOW{1,2}(sw_model_interped);
    pc_sw_obs = model.fluid.pcOW{1,2}(sw_obs_interped);
    

    if plots
        figure('Name', 'F Factors calculations','NumberTitle', 'off', 'WindowStyle','docked');
        tiledlayout('flow')

        nexttile
        p1 = plot(x / Convert('cm'), pc_model / Convert('bar'), 'DisplayName', 'Pc Model');
        title('Homogeneous and heterogeneous case pc');
        hold on
        size_pc = size(pc_model);
        for i = 1:size_pc(1)
            stairs(x / Convert('cm'),pc_sw_obs(i,:)  / Convert('bar'),...
                'Color', p1(i).Color, 'DisplayName', 'Pc Experiment')
        end
        xlabel('Location (cm)')
        ylabel('pressure (bar)')
        grid off
        for i = 1:length(t_obs)
            p1(i).DataTipTemplate.DataTipRows(1).Label = 'Position (cm)';
            p1(i).DataTipTemplate.DataTipRows(2).Label = 'Water Saturation';
            row = dataTipTextRow('Time (hour)', repmat(t(i)...
                / Convert('hr'),1,length(pc_sw_obs)));
            p1(i).DataTipTemplate.DataTipRows(3) = row;
            index = dataTipTextRow('Index', repmat(i...
                ,1,length(pc_sw_obs)));
            p1(i).DataTipTemplate.DataTipRows(4) = index;
        end

        nexttile
        p2 = plot(x / Convert('cm'), pc_model(index_mask, :) / Convert('bar'));
        title('Homogeneous and heterogeneous case pc');
        hold on
        for i = 1:length(index_mask)   
            stairs(x / Convert('cm'), pc_sw_obs(index_mask(i),:) / Convert('bar'),...
                'Color', p2(i).Color)
        end
        xlabel('Location (cm)')
        ylabel('pressure (bar)')
        grid off
        for i = 1:length(index_mask)
            p2(i).DataTipTemplate.DataTipRows(1).Label = 'Position (cm)';
            p2(i).DataTipTemplate.DataTipRows(2).Label = 'Water Saturation';
            row = dataTipTextRow('Time (hour)', repmat(t(i)...
                / Convert('hr'),1,length(pc_sw_obs)));
            p2(i).DataTipTemplate.DataTipRows(3) = row;
            index = dataTipTextRow('Index', repmat(i...
                ,1,length(pc_sw_obs)));
            p2(i).DataTipTemplate.DataTipRows(4) = index;
        end
    end
    
    constant_porosity = model.experiment.rock.poro.value;
    constant_perm = model.experiment.rock.perm.value;
    if isfield(model.experiment.rock.poro, 'porosity_profile')
        porosity_profile_obs = model.experiment.rock.poro.porosity_profile;
        length_obs = porosity_profile_obs(:,1); 
        porosity_obs = porosity_profile_obs(:,2);
        poro_array = interp1(length_obs * Convert('cm'), porosity_obs, x); 
    else
        poro_array = repmat(constant_porosity, length(x), 1);
    end
    
% initialize f with a guess coming from the first pressure profile
    f0 = pc_sw_obs(index_mask(1),:) ./ pc_model(index_mask(1),:);
    obj_fun = @(x) objective_function(x, pc_model(index_mask,:), ...
        pc_sw_obs(index_mask,:), sw_model_interped(index_mask,:),...
        sw_obs_interped(index_mask,:), alpha);
    inequalities = @(iterator) mycon(iterator, poro_array(:), constant_porosity, constant_perm);
    options_fmincon = optimoptions('fmincon','Display','off');
    lb = zeros(length(f0),1); ub = ones(length(f0),1)*2;
    f_factor = fmincon(obj_fun,f0(:),[],[],[],[],lb,ub,inequalities,options_fmincon);
    
% in case one wants to test other optimizers
%     nvars = length(f0);
%     lb = ones(nvars, 1)*0.5;
%     ub = ones(nvars, 1)*1.5;
%     options_swarm = optimoptions('particleswarm','Display','iter');
%     f_factor = particleswarm(obj_fun,nvars,lb,ub,options_swarm);

    if plots
        nexttile
        pc_model_p = plot(model.satfun.sw_pc, model.satfun.pc / Convert('bar'),'r-', 'DisplayName', 'Pc model');
        title('Capillary pressure scaling');
        hold on
        scaled_pc =  f_factor(:)' .* pc_model;
        scale_data_p = plot(sw_obs_interped(index_mask, :), scaled_pc(index_mask, :) / Convert('bar'), 'bs', 'DisplayName', 'Scaled data');
        unscaled_data_p = plot(sw_obs_interped(index_mask, :), pc_model(index_mask, :) / Convert('bar'), 'ks', 'DisplayName', 'Unscaled data');
        legend([pc_model_p(1) scale_data_p(1) unscaled_data_p(1)],{'Pc model', 'Scaled data', 'Unscaled data'})
        xlabel('Water Saturation'); ylabel('Capillary pressure (bar)')
        xlim([min(sw_obs_interped,[],'all') max(sw_obs_interped,[],'all')])

% check with syntetic f_factor
        nexttile
        title('Compare with the synthetic f factor');
        ylim([0.5 1.5]); hold on; legend show; grid on 
        stairs(x / Convert('cm'), f_factor(:), 'DisplayName', 'Optimization results'); 
        xlabel('Location (cm)'); ylabel('Calculated F Factors')
%         rng(3)
%         f_origin = 0.9 + (1.1 - 0.9) .* rand(model.grid.G.cells.num, 1);
%         f_origin(1:2) = 1; f_origin(end-1:end) = 1;
%         stairs(1:model.grid.G.cells.num, f_origin, 'DisplayName', 'Original f factor');
%         title('rms error = ' + string(rms(f_factor(:)-f_origin)))
    end

end

function error = objective_function(f, pc_model, pc_sw_obs, sw_model, sw_obs, alpha)
    f = f(:);
    % alpha defines the weighting trade of between the pressure and the
    % saturation
    difference_p = alpha * (pc_sw_obs - pc_model .* f') ./  pc_sw_obs;
    difference_sw = (1-alpha) * (sw_obs - sw_model) ./  sw_obs;
    difference_total = difference_p + difference_sw;
    difference_total(isnan(difference_total)|isinf(difference_total))=[];
    error = rms(difference_total, 'all');
end
function [c,ceq] = mycon(iterator, poro_array, constant_porosity, constant_perm)
    iterator = iterator(:);
    c(1) = abs(mean(iterator) - 1) - 0.01;
    perm_array = constant_perm .* iterator .^2 .* poro_array ./ constant_porosity;
    k_hm = harmmean(perm_array);
    c(2) = abs(k_hm - constant_perm) ./ constant_perm - 0.1;
    ceq =  [];
end