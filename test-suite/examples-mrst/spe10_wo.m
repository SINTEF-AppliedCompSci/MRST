function setup = spe10_wo(varargin)
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
    % One-line description
    description = ...
        ['SPE10 Model 2 (Christie & Blunt, 2008, SPE Res. Eval. & Eng., ', ...
         'doi: 10.2118/72469-pa) with support for subset of layers'      ];
    % Optional input arguments
    options = struct('layers', []); % Layer subset. Either subset of 1:85, 
                                    % 'tarbert', or 'upper_ness'
    [options, optOnly, setup] = processTestCaseInput(mfilename, options, description, varargin{:});
    % Pick layer subset
    if isempty(options.layers)
        options.layers = 1:85;
    elseif ischar(options.layers)
        switch options.layers
            case 'tarbert'
                options.layers = 1:35;
            case 'upper_ness'
                options.layers = 36:85;
            otherwise
                error(['Optional input argument `layers` must be '  , ...
                       'either a vetor of layer numbers, or one of ', ...
                       'the aliases `tarbert` (= layer 1:35) or '   , ...
                       '`upper_ness` (= layer 35:85)']              );
        end
    else
        assert(min(options.layers) >= 1 && max(options.layers) <= 85);
    end
    if optOnly, setup.options = options; return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil spe10
    % Get state0, model and schedule
    gravity reset on  
    [state0, model, schedule] = setupSPE10_AD(options.extra{:}, 'layers', options.layers);
    % Plotting
    if model.G.cartDims(3) > 2
        view = [120,26];
        proj = 'perspective';
        size = [800, 400];
    else
        view = [0,90];
        proj = 'orthographic';
        size = [500, 800];
    end
    zpba = 0.3.*numel(options.layers)/85;
    plotOptions = {'PlotBoxAspectRatio', [1,1.83,zpba], ...
                   'Projection'        , proj        , ...
                   'View'              , view        , ...
                   'Size'              , size        };
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end