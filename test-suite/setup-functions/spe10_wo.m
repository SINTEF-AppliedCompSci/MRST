function setup = spe10_wo(varargin)
% Setup function for Model 2 from the 10th SPE Comparative Solution Project
%
% SYNOPSIS:
%   setup = spe10_wo('pn1', pv1, ...)
%   setup = spe10_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup function for Model 2 from the 10th SPE Comparative Solution
%   Project (Christie & Blunt, SPE Res. Eval. Eng., 2008, doi:
%   10.2118/72469-PA). This model was initially posed as a benchmark to
%   study upscaling techniques but has widely utilized as an example of a
%   highly heterogeneous reservoir rock.
%
%   The 220x60x85 Cartesian grid in the original model represents part of
%   the Brent group from the North Sea: the 35 top layers are from the
%   shallow-marine Tarbert formation (letter 't' in Brent), whereas the 50
%   bottom layers are from the Upper Ness (letter 'n' in Brent). The
%   permeability is highly heterogeneous with values spanning eight orders
%   of magnitude. 
%
%   The reservoir is produced from a classic five-spot pattern, with one
%   injector in each corner and a producer in the center of the model. The
%   fluid model is a two-phase dead-oil model with weak compressibilities,
%   an oil-water viscosity ratio of approximately 10:1, and quadratic
%   relative permeabilities. It is common to consider the model as
%   incompressible.
%
%   The function can setup subsets of the model consisting of any selection
%   of layers, given by the input parameter:
%      'layers' - a single number or a vector of integer values from the
%                 range 1:85. Alternatively, the strings 'tarbert' and
%                 'upper_ness' will give the respective formation.
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    [options, fullSetup, setup] = processTestCaseInput(mfilename, options, description, varargin{:});
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
    if ~fullSetup, setup.options = options; return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil spe10

    % Get state0, model, and schedule
    gravity reset on  
    [state0, model, schedule] = setupSPE10_AD(options.extra{:}, ...
        'layers', options.layers);

    % Plotting
    if model.G.griddim == 3 && model.G.cartDims(3) > 1
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