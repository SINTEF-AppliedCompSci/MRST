function plot_input_range(model) 
%
% DESCRIPTION: plots the boundaries used for the history matching, this is
%              useful in case of using a saturation function
%              parametrization like LET it is not easy to detect which 
%              combination of parameters corresponds to the min or max
%              kr or pc
%
% SYNOPSIS:
%   plot_input_range(model) 
%
% PARAMETERS:
%   model - struct containing following fields:
%   - history_match: boundaries saved from the settings file and configure
%   module
%
% RETURNS:
%   plot of the boundaries of kr and pc used for history matching
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
    fprintf('plotting input range...\n')
    lb = model.history_match.lb(:);
    ub = model.history_match.ub(:);
    if model.history_match.pc.status % not done
        figure('Name','Capillary Pressure','NumberTitle','off')
        title('Capillary Pressure')
        hold on
    %     plot(sw_lb_pc,lb_pc*1e-5,'r--')
    %     plot(sw_ub_pc,ub_pc*1e-5, 'k--')
    %     plot(sw_pc_sol,solution_pc*1e-5, 'b-*')
        plot(sw_x0_pc,x0_pc*1e-5,'r-o', 'MarkerSize',2)
        plot(sw_sim_result_pc,sim_result_pc*1e-5,'b-o', 'MarkerSize',2)
        axis tight
    %     ylim([0 1])
        xlim([0 1])
    %     legend('Lower bound', 'Upper bound', 'True Value', 'Initial Point', 'Simulation Result')
        legend('Initial Point', 'Simulation Result')
        ylabel('Capillary Pressure (bar)')
        xlabel('Water Saturation (fraction)')
    end
    
    if model.history_match.kr.status 
        figure('Name','Relative Permeability','NumberTitle','off')
        title('Relative Permeability')
        hold on; ax_kr = gca;        
        ylim([0 1])
        xlim([0 1])
        ylabel('Relative Permeability')
        xlabel('Water Saturation (fraction)')
        nPoints = 3; allVecs = {};
        x = linspace(0,1,nPoints);
        spaced_matrix = lb + x.*(ub - lb);
        spaced_matrix_size = size(spaced_matrix);
        for i=1:spaced_matrix_size(1)
            allVecs{end+1} = spaced_matrix(i,:);
        end
        sub = cell(1,numel(allVecs));
        [sub{:}] = ndgrid(allVecs{:});
        sub = cellfun(@(x)x(:),sub,'UniformOutput', false);
        % allPerms is [m x n] matrix of m permutations of n vectors
        % m should equal prod(cellfun(@numel,allVecs))
        % n should equal numel(allVecs)
        allPerms = unique(cell2mat(sub),'rows');
        allPerms_size = size(allPerms);
        for i = 1:allPerms_size(1)
            [Sw,krw,kro,~,~] = saturation_tables(allPerms(i,:),model);
            plot(ax_kr,Sw,krw, 'b--', Sw, kro, 'g--')
%             if rem(i, 1e4) == 0; pause(1e-6); end
        end
    end
    
    function [Sw,krw,kro,Sw_pc,pc] = saturation_tables(x,model)
        % build the saturation functions
        params_no_kr = get_params_no_kr(model);
        model = Create_pc_history_match(x, model, params_no_kr);
        model = Create_kr_history_match(x, model);
        
        Sw = model.satfun.sw_kr;
        Sw_pc = model.satfun.sw_pc;
        pc = model.satfun.pc;
        krw = model.satfun.krw;
        kro = model.satfun.kro;
    end
end