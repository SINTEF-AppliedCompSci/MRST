function plotReportIterations(report, schedule)
% Plot nonlinear convergence behavior for simulateScheduleAD output
%
% SYNOPSIS:
%   plotReportIterations(report, schedule)
%
% DESCRIPTION:
%   Creates a plot in the current figure of the nonlinear convergence
%   behavior of a report and schedule pair. Timesteps and ministeps are
%   visualized, along with green and red boxes that indicate successful
%   and wasted iterations respectively.
%
% REQUIRED PARAMETERS:
%   report   - Report as returned by `simulateScheduleAD`.
%
%   schedule - The schedule passed to `simulateScheduleAD` to produce the
%              report.
%
% RETURNS:
%   Nothing.


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

    if nargin == 1
        schedule = [];
    end
    nctrl = numel(report.ControlstepReports);
    nsubstep = cellfun(@(x) numel(x.StepReports), report.ControlstepReports);
    substepNo = rldecode((1:nctrl)', nsubstep);
    
    ministeps = cellfun(@(x) vertcat(x.StepReports{:}), ...
                report.ControlstepReports, 'UniformOutput', false);
    ministeps = vertcat(ministeps{:});
    
    ok = vertcat(ministeps.Converged);
    its = vertcat(ministeps.Iterations);
    
    T_mini = vertcat(ministeps.LocalTime);
    T_ctrl = [0; report.ReservoirTime(1:end-1)];
    
    tunit = day;
    % Timesteps of all ministeps (including those who failed)
    T_all = (T_mini + T_ctrl(substepNo))/tunit;
    % The actual timesteps that converged
    T = [0; T_mini(ok) + T_ctrl(substepNo(ok))]/tunit;
    newplot;
    hold on
    
    nsub = numel(T_all);
    steps = (1:nsub)';
    
    % Plot the ministeps
    covered = false(nsub, 1);
    stepSize = zeros(nsub, 1);
    for i = 2:numel(T)
        t0 = T(i-1);
        t1 = T(i);

        substeps = steps(T_all > t0 & T_all <= t1);
        
        notConverged = ~ok(substeps);
        if any(notConverged)
            good = substeps(~notConverged);
            bad = substeps(notConverged);
        else
            assert(numel(substeps) <= 1);
            good = substeps;
            bad = [];
        end
        
        offset = 0;
        if any(good)
            plotBar(t0, t1, 0, its(good), [0, 178, 60]/255);
            offset = its(good);
        end
        
        for j = 1:numel(bad)
            ib = its(bad(j));
            plotBar(t0, t1, offset, offset + ib, [178, 5, 0]/255);
            offset = offset + ib;
            
        end

        covered(substeps) = true;
        stepSize(i) = offset;
    end

    T_sub = [0; report.ReservoirTime]/tunit;
    
    % Plot changes in controls
    if ~isempty(schedule)
        all_ctrl = schedule.step.control;
        ctrl = unique(all_ctrl);
        
        if numel(ctrl) > 1
            colors = jet(max(ctrl));
            c = all_ctrl(1);
            x = 0;
            y = max(stepSize);
            for i = 2:numel(all_ctrl)
                if i == numel(all_ctrl) || all_ctrl(i+1) ~= c
                    xn = T_sub(i+1);
                    plotBar(x, xn, 0, y, [1, 1, 1], 'facecolor', 'none', 'edgecolor', colors(c, :), 'linewidth', 3);
                    x = xn;
                    if i < numel(all_ctrl)
                        c = all_ctrl(i+1);
                    end
                end
            end
        end
    end
    
    
    % Plot control steps
    for i = 1:nctrl
        plotBar(T_sub(i), T_sub(i+1), 0, max(stepSize), [1, 1, 1], ...
            'FaceColor', 'none', 'EdgeColor', [.65, .65, .65], 'LineWidth', 1, 'LineStyle', '-');
    end
    ylabel('Iterations used in ministep');
    xlabel('Time [days]')
end

function h = plotBar(x1, x2, y0, y1, varargin)
    X = [x1, x1, x2, x2];
    Y = [y0, y1, y1, y0];
    h = patch(X, Y, varargin{:});
end