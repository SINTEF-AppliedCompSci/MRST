function [schedule, report, isAltered] = controlLogicFunc(state, schedule, report, cnum, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('wcutLim',                   [], ...
             'rateLim',                   [], ...
             'closeWellOnSignChange', false);

opt = merge_options(opt, varargin{:});
curc = schedule.step.control(cnum);
[W, ws]  = deal(schedule.control(curc).W, state.wellSol);
[nc, nw] = deal(numel(schedule.control), numel(W));
isAltered = false;
for k  =1:nw
    % check watercut in producers
    wclim = opt.wcutLim;
    if ~isempty(wclim) && W(k).sign < 0 && W(k).status > 0
        wc = ws(k).qWs/(ws(k).qWs + ws(k).qOs);
        if isfinite(wc) && wc >= wclim
            isAltered = true;
            fprintf('Well %s is shut down due to high water-cut (%3.2f)\n', W(k).name, wc);
            for c = curc:nc
                schedule.control(c).W(k).status = false;
            end
        end
    end
    % check rates
    ratelim = opt.rateLim;
    if ~isempty(ratelim)
        rate = (ws(k).qWs + ws(k).qOs);
        if abs(rate) < ratelim && W(k).status > 0
            isAltered = true;
            fprintf('Well %s is shut down due to small rate (%3.2f)\n', W(k).name, rate);
            schedule.control(curc).W(k).status = false;
        end
    end
    % check well sign change
    if opt.closeWellOnSignChange
        rate = (ws(k).qWs + ws(k).qOs);
        if (sign(rate) ~= W(k).sign)
            isAltered = true;
            if W(k).sign > 0
                fprintf('Injector %s is shut down due to negative rate: (%3.2f)\n', W(k).name, rate);
            else
                fprintf('Producer %s is shut down due to positive rate: (%3.2f)\n', W(k).name, rate);
            end
            schedule.control(curc).W(k).status = false;
        end
    end
end
end
