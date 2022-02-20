function plotWellRates(W, data, varargin)
% Simple utility for plotting well rates. Used by optimizeTOF.

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

    if size(data, 2) ~= numel(W)
        data = data .';
    end
    if size(data, 1) < 2
        return
    end
    assert(size(data, 2) == numel(W), 'Input data and well number must match!');
    subs = repmat(1:size(data, 1), 2, 1);
    subs = subs(:);

    t = 1:size(data, 1);
    % Normalize
    data = bsxfun(@rdivide, data, sum(data(end, :), 2));

    area(t(subs(2:end-1)), data(subs(1:end-2), :));
    axis tight
    l = arrayfun(@(x) x.name, W, 'UniformOutput', false);
    caxis([.5 size(data, 2) - .5]);
    if nargin > 2 && any(varargin{1})
        % We were given indices of steps where improvement to the objective
        % function happened, plot them as well!
        hold on
        improvement = varargin{1};
        if islogical(varargin{1})
            improvement = find(improvement);
        end
        ni = numel(improvement);
        plot(improvement, zeros(ni, 1), 'd', 'markerfacecolor', 'g', 'markersize', 12)

        l = [l; 'Improvement'];
    end
    legend(l);
end
