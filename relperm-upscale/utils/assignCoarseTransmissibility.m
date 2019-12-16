function [model, schedule] = assignCoarseTransmissibility(p, model, schedule, T_struct, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('timeMethod',    'mean', ...
                 'timeFunction',  [], ...
                 'phaseMethod',   'none', ...
                 'phaseFunction', [], ...
                 'replaceNegatives', 'zero', ...
                 'replaceNans',      'zero', ...
                 'fallbackTransmissibility', []);
    opt = merge_options(opt, varargin{:});
    
    
    r = T_struct.reservoir;
    T_r = processConnection(schedule.step.val, r.T, r.T_ph, opt);
    T_r = replaceNegativesAndNans(T_r, opt);
    
    model = upscaleModelTPFA(model, p, 'transCoarse', T_r, 'neighborship', r.N);
    schedule = upscaleSchedule(model, schedule);
    if nargout > 1
        w = T_struct.wells;
        for ctrl = 1:numel(schedule.control)
            WI = processConnection(schedule.step.val, w.WI, w.WI_ph, opt);
            if ~all(WI>=0)
                warning('Upscaling process produced some negative well transmissibilities which are set to zero');
            end
            WI(WI<0) = 0;
            for wNo = 1:numel(schedule.control(ctrl).W)
                schedule.control(ctrl).W(wNo).WI = WI(w.perf2well == wNo);
            end
        end
    end
end

function T = replaceNegativesAndNans(T, opt)
if ~isempty(opt.fallbackTransmissibility)
    [opt.replaceNegatives, opt.replaceNans] = deal('fallback');
end
switch opt.replaceNegatives
    case 'none'
        % don't do anything
    case 'zero'
        T(T<0) = 0;
    case 'fallback'
        assert(~isempty(opt.fallbackTransmissibility))
        T(T<0) = opt.fallbackTransmissibility(T<0);
end
switch opt.replaceNans
    case 'no'
        % don't do anything
    case 'zero'
        T(isnan(T)) = 0;
    case 'fallback'
        assert(~isempty(opt.fallbackTransmissibility))
        T(isnan(T)) = opt.fallbackTransmissibility(isnan(T));
end
end

% function T = fixTrans(T)
%     T(~isfinite(T)) = 0;
%     if any(T > 0)
%         T(T<0) = 1e-3*min(T(T>0));
%     else
%         T(T<0) = 0;
%         warning('All transmissibilities are negative!');
%     end
% end

function T = processConnection(dt, T, T_ph, opt)
    nph = numel(T_ph);
    
    %T = fixTrans(T);
    %for i = 1:nph
    %    T_ph{i} = fixTrans(T_ph{i});
    %end
    if isempty(opt.phaseFunction)
        switch lower(opt.phaseMethod)
            case 'max'
                T = -inf;
                for i = 1:nph
                    T = max(T, T_ph{i}, 'omitnan');
                end
            case 'min'
                T = inf;
                for i = 1:nph
                    T = min(T, T_ph{i},'omitnan');
                end
            case 'mean'
                mean(horzcat(T_ph{:}), 2, 'omitnan');
                T  = zeros(size(T_ph{1}));
                nn = 0;
                for i = 1:nph
                    ix    = ~isnan(T_ph{i});
                    T(ix) = T(ix) + T_ph{i}(ix);
                    nn = nn + ix;
                end
                T = T./nn;
            case 'none'
                % Just use total trans
            otherwise
                error(['Unknown matching method "', opt.phaseMethod, '".']);
        end
    else
        T = opt.phaseFunction(T, T_ph, dt);
    end
    
    if isempty(opt.timeFunction)
        switch lower(opt.timeMethod)
            case 'max'
                T = max(T, [], 2, 'omitnan');
            case 'min'
                T = min(T, [], 2, 'omitnan');
            case 'mean'
                T = bsxfun(@times, T, dt');
                T = sum(T, 2, 'omitnan')./((~isnan(T))*dt);
            otherwise
                error(['Unknown matching method "', opt.timeMethod, '".']);
        end
    else
        T = opt.timeMethod(T, T_ph, dt);
    end
end
