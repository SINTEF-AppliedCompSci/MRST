function [description, options, state0, model, schedule, plotOptions] = spe1_bo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`,

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
    % One-line description
    description = ...
        ['SPE1 benchmark, see Odeh, A.S., "Comparison of Solutions to '            , ...
         'a Three-Dimensional Black-Oil Reservoir Simulation Problem.", '          , ...
         'J. Pet. Technol. 33 (1): 13–25. SPE-9723-PA. 1981, doi: 10.2118/9723-PA '];
    % Optional input arguments
    options = struct();
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil deckformat
    % Get initial state, model and schedule
    pth = getDatasetPath('spe1');
    fn  = fullfile(pth, 'BENCH_SPE1.DATA');
    [state0, model, schedule] = initEclipseProblemAD(fn);
    % Set plot options
    plotOptions = {'PlotBoxAspectRatio', [1,1,0.25], ...
                   'View'              , [-15, 20] };
end