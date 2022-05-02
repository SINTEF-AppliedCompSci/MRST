function state = setWellValuesFromGroupControls(model, state0, state, dt, drivingForces)
% Set well values for all wells that operates under group control
%

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

    % Check if groups are present
    groups = drivingForces.groups;
    ng     = numel(groups);
    if ng == 0, return; end
    % Compute well potential
    pot = model.computeWellPotential(state);

    groups = processGroups(model, groups, state, pot);
    ng     = numel(groups);
    % Loop through groups and update well controls
    for g = 1:ng
        group = groups{g};
        switch group.type
            case 'none'
                % Group does not have a control - nothing to do
                continue
            case 'rate'
                state = setGroupWellRates(model, group, state, pot);
            case 'temperature'
                if ~model.ReservoirModel.thermal, continue; end
                state = setGroupWellTemperatures(model, group, state);
            otherwise
                error('Group control type not supported');
        end
    end

end

%-----------------------------------------------------------------%
function state = setGroupWellRates(model, group, state, pot)
% Set well rates based on group target

    % Get well solutions
    ws = state.wellSol;
    map = model.getProps(state, 'FacilityWellMapping');
    if ~any(map.active), return; end
    % Get group wellSols
    mask = getGroupMask(model, state, group.name);
    wsg  = ws(mask);
    % Get well rates and temperatures for all wells in group
    q = model.getWellRates(state);
    T = model.getWellTemperatures(state);
    if ischar(group.val) || iscell(group.val)
        % Control value is the name of another group - this means
        % that we aim at a group rate equal to that groups rate,
        % with opposite sign
        maskp = getGroupMask(model, state, group.val);
        qtot  = -sum(q(maskp));
        Ttot  = sum(q(maskp).*T(maskp))./sum(q(maskp));
    elseif isnumeric(group.val)
        % Control value is a numeric variable
        qtot = group.val;
        Ttot = group.T;
    end
    if isfield(group, 'frac')
        qtot = qtot*group.frac;
    end
    % Loop through wells in group
    active = true(numel(wsg),1);
    wix    = find(mask);
    for i = 1:numel(wsg)
        if wsg(i).sign > 0
            % Set target temperature is injector
            state.FacilityFluxProps.FacilityWellMapping.W(wix(i)).T = Ttot;
        end
        if ~strcmpi(wsg(i).type, 'rate')
            % Well is not operating at rate. This should mean that
            % it has reached its capacity, so we extract its rate
            % from the total goal
            qw   = q(wix(i));
            qtot = qtot - qw;
            active(i) = false;
        end
    end
    % Distribute total target group rate among all wells not
    % working at maximum capacity
    wsg_act = wsg(active);
    pot     = pot(mask); pot = pot(active);
    w       = pot./sum(pot);
    for i = 1:numel(wsg_act)
        wsg_act(i).val = qtot.*w(i);
    end
    % Set updated wellSol filed to state
    wsg(active)   = wsg_act;
    ws(mask)      = wsg;
    state.wellSol = ws;

end

%-----------------------------------------------------------------%
function [state, Tw] = setGroupWellTemperatures(model, group, state)
% Set well rates based on group target

    % Get well solutions
    ws = state.wellSol;
    map = model.getProps(state, 'FacilityWellMapping');
    if ~any(map.active), return; end
    % Get group wellSols
    mask = getGroupMask(model, state, group.name);
    wsg  = ws(mask);
    % Get well rates and temperatures for all wells in group
    q = model.getWellRates(state);
    T = model.getWellTemperatures(state);
    if ischar(group.val) || iscell(group.val)
        % Control value is the name of another group - this means
        % that we aim at a group rate equal to that groups temperature,
        maskp = getGroupMask(model, state, group.val);
        Ttot  = sum(q(maskp).*T(maskp))./sum(q(maskp));
        if isfield(group, 'ctrlFn')
            Ttot = group.ctrlFn(Ttot);
        end
    elseif isnumeric(group.val)
        % Control value is a numeric variable
        Ttot = group.T;
    end
    % Loop through wells in group
    wix    = find(mask);
    Tw     = model.AutoDiffBackend.convertToAD(vertcat(ws.T), q);
    for i = 1:numel(wsg)
        if wsg(i).sign > 0
            % Set target temperature is injector
            state.FacilityFluxProps.FacilityWellMapping.W(wix(i)).T = Ttot;
        end
    end

end