function [description, options, state0, model, schedule, plotOptions] = egg_wo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

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

    % Description
    description = ['Egg model, see Jansen, J. D., et al. '       , ...
                   '"The egg model–a geological ensemble for '   , ...
                   'reservoir simulation." '                     , ...
                   'Geoscience Data Journal 1.2 (2014): 192-195.'];
    % Optional input arguments
    options = struct('realization', 0); % Realization [0, 100]
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Get deck
    deck = getDeckEGG('realization', options.realization);
    % Initialize MRST problem from deck
    [state0, model, schedule] = initEclipseProblemAD(deck);
    % Plotting
    plotOptions = {'View'              , [-30, 20] , ...
                   'PlotBoxAspectRatio', [1,1,0.2] , ...
                   'Box'               , true      , ...
                   'Size'              , [800, 500]};
end
