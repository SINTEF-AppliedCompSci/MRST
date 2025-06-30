function [p_opt, history] = improvedCalibrationGame(pvec_1,pvec_2, objh2, objh2_new, varargin)
    % Parameters with defaults
    maxIt_alt = 5;
    maxIt = 5;
    objTol = 1e-6;
    paramTol = 1e-4;
    stallIter = 3;
    scale_init = 0.7;
    
    % Initialize
    p_opt = pvec_1;
    history_main = struct('val', [], 'u', {{}}, 'lambda', [], 'rho', [], 'nIt', [], 'pg', [], 'du', []);
    no_improvement = 0;    
    for i = 1:maxIt_alt
        % Adaptive scaling based on rho (gain ratio)
        if i > 1 && ~isnan(history_main.rho(end))
            current_scale = scale_init * min(1.5, max(0.1, history_main.rho(end)));
        else
            current_scale = scale_init;
        end
        
        % First optimization step with adaptive LM parameters
        if i > 1 && length(history_main.lambda) >= 1
            lambda_init = history_main.lambda(end);
        else
            lambda_init = 0.01;
        end
        
        [~, p_new, hist_new] = unitBoxLM(p_opt, objh2_new, ...
            'maxIt', maxIt, 'lambda0', lambda_init);
        
        % Adaptive mixing using gradient information if available
        if isfield(hist_new, 'pg') && ~isempty(hist_new.pg)
            grad_norm = hist_new.pg(end);
            mix_factor = min(0.3, 0.05 + 0.25/(1 + grad_norm));
            p_new = mix_factor.*p_new + (1-mix_factor).*p_opt;
        else
            p_new = current_scale.*p_new + (1-current_scale).*p_opt;
        end
        
        history_main = appendHistory(history_main, hist_new);
        
        % Second optimization step with direction-aware mixing
        [~, p_new2, hist2] = unitBoxLM(p_new, objh2, 'maxIt', maxIt);
        
        % Direction-aware mixing using du (step direction) if available
        if isfield(hist2, 'du') && ~isempty(hist2.du) && ~isnan(hist2.du(end))
            direction_factor = min(1, max(0, hist2.du(end)));
            p_new2 = (current_scale*direction_factor).*p_new2 + ...
                    (1-current_scale*direction_factor).*p_new;
        else
            p_new2 = current_scale.*p_new2 + (1-current_scale).*p_new;
        end
        
        history_main = appendHistory(history_main, hist2);
        
        % Enhanced convergence checking using multiple criteria
        [converged, no_improvement] = checkEnhancedConvergence(...
            history_main, p_opt, p_new2, objTol, paramTol, no_improvement, stallIter);
        
        if converged
            break;
        end
        
        % Update for next iteration
        p_opt = p_new2;
        
        % Adaptive maxIt reduction if making good progress
        if i > 2 && mean(diff(history_main.val(end-3:end))) < -objTol/10
            maxIt = max(3, ceil(maxIt * 0.8));
        end
    end
    
    p_opt = p_new2;
    history = history_main;
end

function [converged, no_improvement] = checkEnhancedConvergence(...
        history, p_old, p_new, objTol, paramTol, no_improvement, stallIter)
    
    if length(history.val) < 2
        converged = false;
        return;
    end
    
    % Calculate multiple convergence metrics
    obj_change = abs(history.val(end) - history.val(end-1));
    rel_obj_change = obj_change / abs(history.val(end-1));
    param_change = norm(p_new - p_old)/norm(p_old);
    
    % Additional metrics from history
    if isfield(history, 'pg') && length(history.pg) >= 1
        grad_norm = history.pg(end);
    else
        grad_norm = inf;
    end
    
    if isfield(history, 'rho') && length(history.rho) >= 1 && ~isnan(history.rho(end))
        rho_value = history.rho(end);
    else
        rho_value = 1;
    end
    
    % Composite convergence criteria
    criterion1 = (rel_obj_change < objTol) && (param_change < paramTol);
    criterion2 = (grad_norm < 1e-3);
    criterion3 = (rho_value < 0.1) && (obj_change < objTol);
    
    if criterion1 || criterion2 || criterion3
        no_improvement = no_improvement + 1;
        converged = (no_improvement >= stallIter);
    else
        no_improvement = 0;
        converged = false;
    end
end

function history = appendHistory(history, new_history)
    % Append all fields from new_history to history
    fields = fieldnames(new_history);
    for k = 1:numel(fields)
        field = fields{k};
        if isfield(history, field)
            if iscell(history.(field))
                history.(field) = [history.(field) new_history.(field)];
            else
                history.(field) = [history.(field), new_history.(field)];
            end
        else
            history.(field) = new_history.(field);
        end
    end
end