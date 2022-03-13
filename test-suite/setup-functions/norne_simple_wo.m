function setup = norne_simple_wo(varargin)
% Setup function for a water-oil ensemble version of the Norne field
%
% SYNOPSIS:
%   setup = norne_simple_wo('pn1', pv1, ...)
%   setup = norne_simple_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a simplified water-oil ensemble model of the Norne field. The
%   reservoir geometry is from the full simulation model published by the
%   Open Porous Media (OPM) project. Petrophysical data are generated using
%   geostatistical routines developed in a paper by Lorentzen et al, SPE J,
%   2019, doi: 10.2118/194205-PA, but without inclusion of some of the
%   random parameters that modify connectivity across faults.
%
%   The production scenario is greatly simplified compared with the real
%   field and describes immiscible waterflooding into a reservoir that
%   initially is fully filled with oil. On purpose, the wells are placed
%   and controlled so that one of the producers will dominate the overall
%   field production.
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial, norne_field_bo

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
    description = ['The Norne field model used in the paper: '  , ...
                   'Flow Diagnostics for Model Ensembles, '     , ...
                   'https://doi.org/10.3997/2214-4609.202035133'];

    % Optional input arguments
    options = struct('realization', 1);
    [options, fullSetup, setup] = processTestCaseInput(mfilename, options, description, varargin{:});    
    if ~fullSetup, return; end
    
    % Module dependencies
    require ad-core ad-props ad-blackoil
    
    % Get ensemble
    out      = setupNorneFn(options.realization); 
    state0   = out.state0;
    model    = out.model;
    fluid = initSimpleADIFluid('phases', 'WO'                        , ...
                               'mu'    , [   1,   5]*centi*poise     , ...
                               'rho'   , [1014, 859]*kilogram/meter^3, ...
                               'n'     , [   2,   2]                 );
    model = GenericBlackOilModel(model.G, model.rock, fluid, ...
                                 'gas'      , false          , ...
                                 'inputdata', model.inputdata);
    schedule = out.schedule;

    % Plotting
    plotOptions = {'View'              , [85, 30]  , ...
                   'PlotBoxAspectRatio', [1,1,0.4] , ...
                   'Size'              , [700, 600]};
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end