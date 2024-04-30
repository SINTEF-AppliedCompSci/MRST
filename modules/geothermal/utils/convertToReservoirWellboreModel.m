function setup = convertToReservoirWellboreModel(setup, varargin)
% Convert a setup with geothermal model to composite model with WellboreModel for wells
%
% SYNOPSIS:
%   setup = convertToReservoirWellbore(setup)
%   setup = convertToReservoirWellbore(setup, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   setup  - Test case setup. Can either be an instance of `TestCase`, or a
%            struct with at least the fields `model` (reservoir model),
%            `schedule` (simulation schedule), and `state0` (initial
%            state).
%  
% RETURNS:
%   setup - Updated setup with composite model where the wells in the
%           schedule have been converted to an instance of WellboreModel.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    % Optional input arguments
    opt = struct('name', []);
    [opt, varargin] = merge_options(opt, varargin{:});
    if isempty(opt.name), opt.name = [setup.name, '_wb']; end
    
    % Get reservoir model and construct wellbore model
    reservoirModel = setup.model;
    W = setup.schedule.control(1).W;    
    wellboreModel = WellboreModel(reservoirModel, W, varargin{:});
    % Make composite model from the two
    model = CompositeModel({reservoirModel, wellboreModel}, ...
                              'names', {'Reservoir', 'Wellbore'});
                          
    % Define couplings between wellbore and 
    % Mass flux
    coupling = WellboreReservoirComponentFlux(model, ...
        'Reservoir', 'Wellbore');
    model = model.setCouplingTerm(coupling);
    if reservoirModel.thermal
        % Heat flux
        coupling = WellboreReservoirHeatFlux(model, ...
        'Reservoir', 'Wellbore');
        model = model.setCouplingTerm(coupling);
    end
    
    % Set up initial state
    state0 = struct();
    if isfield(setup.state0, 'wellSol')
        setup.state0 = rmfield(setup.state0, 'wellSol');
    end
    [state0.Reservoir, state0.Wellbore] = deal(setup.state0);
    wc = wellboreModel.G.cells.global;
    state0.Wellbore.pressure = state0.Wellbore.pressure(wc);
    state0.Wellbore.s        = state0.Wellbore.s(wc);
    state0.Wellbore.T        = state0.Wellbore.T(wc);
    
    % Convert schedule
    schedule = setup.schedule;
    control = schedule.control;
    for i = 1:numel(control)
        [rctrl, wctrl] = deal(control(i));
        [rctrl.W, rctrl.groups] = deal([]);
        if ~isfield(rctrl, 'bc' ), rctrl.bc  = []; end
        if ~isfield(rctrl, 'src'), rctrl.src = []; end
        [wctrl.bc, wctrl.src] = deal([]);
        schedule.control(i).Reservoir = rctrl;
        schedule.control(i).Wellbore = wctrl;
    end
    if isfield(schedule.control, 'W')
        schedule.control = rmfield(schedule.control, 'W');
    end
    if isfield(schedule.control, 'src')
        schedule.control = rmfield(schedule.control, 'src');
    end
    if isfield(schedule.control, 'bc')
        schedule.control = rmfield(schedule.control, 'bc');
    end
    if isfield(schedule.control, 'groups')
        schedule.control = rmfield(schedule.control, 'groups');
    end
    
    % Update setup
    setup.state0   = state0;
    setup.model    = model;
    setup.schedule = schedule;
    setup.name     = opt.name;
    
    % Set visualization grid
    setup.visualizationGrid = model.submodels.Reservoir.G;
    
end