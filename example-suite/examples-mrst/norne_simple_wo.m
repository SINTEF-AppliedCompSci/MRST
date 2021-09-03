function  [description, options, state0, model, schedule, plotOptions] = norne_simple_wo(varargin)
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
    description ...
            = ['The Norne field model used in the paper:  ', ...
               'Flow Diagnostics for Model Ensembles , https://doi.org/10.3997/2214-4609.202035133'];
    % Optional input arguments
    options = struct('realization', 1);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
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
end