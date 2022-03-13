function setup = spe9_bo(varargin)
% Setup function for the ninth SPE Comparative Solution Project (SPE 9)
%
% SYNOPSIS:
%   setup = spe9_bo('pn1', pv1, ...)
%   setup = spe9_bo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup for the ninth SPE Comparative Solution Project (J.E. Killough,
%   Ninth SPE Comparative Solution Project: A reexamination of black-oil
%   simulation, 1995 SPE Reservoir Simulation Symposium, doi:
%   10.2118/29110-MS).
%
%   The test case models a water-injection problem in a highly
%   heterogeneous reservoir consisting of fourteen unique layers having
%   isotropic, lognormally distributed permeability spanning six orders of
%   magnitude. The reservoir geometry is described by a 24x25x15 regular
%   grid with a 10-degree dip in the x-direction. The live-gas, black-oil
%   model has several challenging features including strong kinks and steep
%   gradients in the relative permeability and the water-oil capillary
%   pressure curves.  Other challenges include transition from
%   undersaturated to saturated state, abrupt changes in prescribed
%   production rates, and dynamic switching between rate and bhp control
%   for the producers.
%
%   The model is so large that we replace the standard AD backend in MRST
%   by a more optimized row-major version.
%
%   Simulation of this case is discussed in more detail in section 12.4.4
%   of the first MRST textbook (Lie, Cambridge University Press, 2019, doi:
%   10.1017/9781108591416).
%
%   Configurable parameters: none.
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial, spe1_bo, spe3_bo.


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
        ['SPE9 benchmark, see Killough, J.E., "Ninth SPE comparative  ' , ...
         'solution project: A reexamination of black-oil simulation.", ', ...
         'SPE Reservoir Simulation Symposium, 1995, SPE 29110-MS, '     , ...
         'doi: 10.2118/29110-MS'                                        ];

    % Optional input arguments
    options = struct();
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil deckformat
    
    % Get initial state, model, and schedule
    pth = getDatasetPath('spe9');
    fn  = fullfile(pth, 'BENCH_SPE9.DATA');
    [state0, model, schedule] = ...
        initEclipseProblemAD(fn, 'useMex', true, 'rowMajorAD', true);

    % Set plot options
    plotOptions = {'PlotBoxAspectRatio', [1,1,0.25], ...
                   'View'              , [120, 20] };
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end