function setup = egg_wo(varargin)
% Setup function for the Egg benchmark models
%
% SYNOPSIS:
%   setup = egg_wo('pn1', pv1, ...)
%   setup = egg_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a single model realization from the Egg ensemble; see J.D.
%   Jansen et al., The egg model - a geological ensemble for reservoir
%   simulation, Geoscience Data Journal 1.2 (2014): 192-195. The Egg
%   ensemble contains 101 realizations of a small fluivial reservoir with a
%   perimeter the loosely resembles an egg (hence the name). Each
%   realization  is represented on a 60x60x7 corner-point grid, out of
%   which 18 553 cells are active.
%   
%   The reservoir is produced under waterflooding conditions with eight
%   water injectors operating and constant rate and four producers
%   operating at constant bottom-hole pressure. The oil-water system
%   follows a relatively simple black-oil model with weak compressibility,
%   cubic and quartic relative permeabilities curves, and 5 cP and 1 cP
%   viscosities for the oil and water phases, respectively.
%
%   The only configurable parameter is:
%      'realization' - realization number (default: 0, range: 0 to 100)
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
    % Description
    description = ['Egg model, see Jansen, J. D., et al. '       , ...
                   '"The egg model - a geological ensemble for '   , ...
                   'reservoir simulation." '                     , ...
                   'Geoscience Data Journal 1.2 (2014): 192-195.'];
    % Optional input arguments
    options = struct('realization', 0); % Realization [0, 100]
    
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end

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

    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end
