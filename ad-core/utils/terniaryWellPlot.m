function varargout = terniaryWellPlot(wellSols, T, ix, varargin)
% Plot well curves (water, gas, oil and optionally BHP) for wellSols
%
% SYNOPSIS:
%   % Plot well #3, with timesteps on xaxis
%   terniaryWellPlot(wellSols, time, 3);
%   % Plot all wells
%   terniaryWellPlot(wellSols)
%
% DESCRIPTION:
%   This function is tailored towards three-phase simulation and is capable
%   of producing plots that include production rates for each well for each
%   phase (water, oil gas) and optionally also bottom hole pressures as a
%   separate axis. One figure is produced per well requested.
%
% PARAMETERS:
%   wellSols - Cell array of NSTEP by 1, each containing a uniform struct
%              array of well solution structures. For example, the first
%              output from simulateScheduleAD. Can also be a cell array of
%              such cell arrays, for comparing multiple simulation
%              scenarios.
%
%   time     - (OPTIONAL) The time for each timestep. If not provided, the
%              plotter will use step number as the x axis intead. If
%              wellSols is a cell array of multiple datasets, time should
%              also be a cell array, provided not all datasets use the
%              same timesteps.
%
%   ix       - (OPTIONAL) A list of indices to plot, or a single string
%              corresponding to the name of a specific well. The default is
%              all wells.
%
% OPTIONAL PARAMETERS:
%   'plotBHP'  - Boolean indicating if BHP is to be plotted. Defaults to
%                enabled.
% 
% RETURNS:
%   fh         - Figure handles to all figures that were created.
%
% SEE ALSO:
%   `simulateScheduleAD`, `plotWellSols`

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

    opt = struct('plotBHP', true);
    opt = merge_options(opt, varargin{:});
    
    ws0 = wellSols{end};
    if nargin < 2
        T = 1:numel(wellSols);
        xname = 'Step #';
    else
        T = T/year;
        xname = 'Time [year]';
    end
    if nargin < 3
        ix = 1:numel(ws0);
    end
    
    if ischar(ix)
        ix = find(strcmpi({ws0.name}, ix));
    end
    
    flds = {'bhp', 'qWs', 'qOs', 'qGs'};
    [data, names] = getWellOutput(wellSols, flds, ix);
    figs = zeros(numel(ix), 1);
    for i = 1:numel(ix)
        figs(i) = figure('Color', 'w');
        hold on
        qw = squeeze(data(:, i, 2:4));
        
        c = 'brg';
        fac = 1;
        for j = 1:3
            if j == 3, fac = 1000; end
            d = abs(qw(:, j))*day;
            plot(T, d/fac, [c(j), '-'], 'linewidth', 2)
        end
        ax = gca;
        set(ax, 'Color','none');
        xlabel(xname);
        legend('Water', 'Oil', 'Gas', 'Location', 'NorthWest');
        ylabel('Rates [m^3/day] (liquid),  [1000m^3/day] (gas)');
        if opt.plotBHP
            axpos = get(ax, 'Position');
            ax2 = axes('Position',axpos,...
                        'XAxisLocation','top',...
                        'YAxisLocation','right',...
                        'Color','none');
            hold on
            plot(T, data(:, i, 1)/mega, 'k-', 'Parent', ax2, 'linewidth', 2);
            set(ax2, 'XTickLabel', []);
            ylabel('Bottom hole pressure [MPa]');
            legend('BHP', 'Location', 'NorthEast');
            axes(ax);
        end
        grid on
        
        
        % Title - change depending on well sign
        sgn = ws0(i).sign;
        if sgn == -1
            title(['Production for ', names{i}]);
        elseif sgn == 1
            title(['Injection for ', names{i}]);
        else
            % Who knows what this well is? Hopefully the user does.
            title(['Well curves for ', names{i}]);
        end
    end
    if nargout
        varargout{1} = figs;
    end
end
