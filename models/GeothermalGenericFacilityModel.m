classdef GeothermalGenericFacilityModel < GenericFacilityModel
%Generic facility model for geothermal simulations
    
    properties
        thermalFormulation   = 'temperature';
        implicitConnectionDP = false;
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, useDefaults)
        % Constructor
        
            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end
            % Set up state function groupings using pareng, and add
            % model-specific functionality
            model = setupStateFunctionGroupings@GenericFacilityModel(model, useDefaults);
            ffd = model.FacilityFlowDiscretization;
            ffd = ffd.setStateFunction('AdvectiveHeatFlux' , WellAdvectiveHeatFlux(model) );
            ffd = ffd.setStateFunction('ConductiveHeatFlux', WellConductiveHeatFlux(model));
            ffd = ffd.setStateFunction('HeatFlux'          , HeatFlux(model)              );
            model.FacilityFlowDiscretization = ffd;
            
        end
        
        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
        % Initialize AD state from double state
            
            % Parent model handles standard AD initialization
            state = initStateAD@GenericFacilityModel(model, state, vars, names, origin);
            % If connection pressure drops are  explicit, we are done
            if ~model.implicitConnectionDP, return; end
            % Compute connection pressure drops with derivatives, assuming
            % the wellbore fluid is at uniform pressure (BHP) and
            % temperature (Tw)
            assert(model.ReservoirModel.getNumberOfPhases() == 1,   ...
                ['Implicit connection pressure drop is currently ', ...
                 'only supported for single-phase flow']          );
            % Facility well mapping and gravity
            map    = model.getProp(state, 'FacilityWellMapping');
            g      = norm(model.ReservoirModel.gravity);
            % Get bottom-hole pressure
            isBHP  = strcmpi(state.FacilityState.names, 'bhp');
            bhp    = state.FacilityState.primaryVariables{isBHP};
            % Get well temperature
            isTemp = strcmpi(state.FacilityState.names, 'well_temperature');
            Tw     = state.FacilityState.primaryVariables{isTemp};
            % Get NaCl concentraion if present
            xw     = vertcat(map.W.components);
            Xw     = model.ReservoirModel.getMassFraction(xw);
            cnames = model.ReservoirModel.getComponentNames();
            isNaCl = strcmpi(cnames, 'NaCl');
            if any(isNaCl), Xw = Xw(:,isNaCl); else, Xw = []; end
            % Compute density
            rhow = model.ReservoirModel.fluid.rhoW(bhp, Tw, Xw);
            rhow = rhow(map.perf2well);
            % Compute connection pressure drop for each well
            for i = 1:numel(map.W)
                wellSol = state.wellSol(map.active(i));
                rhowi = rhow(map.perf2well == i);
                wellSol.cdp = rhowi.*g.*map.W(i).dZ;
                state.wellSol(map.active(i)) = wellSol;
            end

        end
        
        %-----------------------------------------------------------------%
        function names = getBasicPrimaryVariableNames(model)
        % Get names of primary variables
        
            % Get parent class primary variable names
            names = getBasicPrimaryVariableNames@GenericFacilityModel(model);
            % Check if we have thermal primary variables
            if strcmpi(model.thermalFormulation, 'none') || ...
               strcmpi(model.primaryVariableSet, 'none')
                return
            end
            % Add well temperature to list of variables
            if model.ReservoirModel.thermal
                names = [names, 'well_temperature'];
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
        % Get primary variables
        
            % Parent model handles almost everything
            [variables, names, map] ...
                = getBasicPrimaryVariables@GenericFacilityModel(model, wellSol);
            if model.ReservoirModel.thermal
                % Add mask indicating where we find the temperature
                map.isTemp = strcmpi(names, 'well_temperature');
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
        % Get value from state based on name
        
            index = ':';
            switch lower(name)
                case 'well_temperature'
                    fn = 'T';
                otherwise
                    [fn, index] = getVariableField@GenericFacilityModel(model, name, varargin{:});
            end
            
        end
        
        %-----------------------------------------------------------------%
        function src = getEnergySources(facility, state)
        % Get thermal energy source terms that goes into the residual
        % equation for conservation of energy
        
            map = facility.getProps(state, 'FacilityWellMapping');
            if isempty(map.W)
                val = [];
            else
                val = facility.getProps(state, 'HeatFlux');
            end
            src = struct('value', {{val}}, 'cells', map.cells);
            
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
        % Get well equations
        
            % If wells are under group control, we must update
            % state.wellSol so that we use correct values in the equations
            state = model.setWellValuesFromGroupControls(state0, state, dt, drivingForces);
            % Call parent model to get mass conservation and control eqs
            [eqs, names, types, state] = getModelEquations@GenericFacilityModel(model, state0, state, dt, drivingForces);
            isTemp = strcmpi(state.FacilityState.names, 'well_temperature');
            if strcmpi(model.thermalFormulation, 'none') || ~any(isTemp) || isempty(eqs)
                return
            end
            % Add control equation for the well temperature
            T   = model.getWellTemperature(state);
            Tw  = state.FacilityState.primaryVariables{isTemp};
            Teq = (Tw - T)./max(norm(value(T),inf), norm(value(Tw),inf));
            eqs = [eqs, {Teq}];
            names = [names, {'temperatureWells'}];
            types = [types, {'perf'}];
            
        end
        
        %-----------------------------------------------------------------%
        function T = getWellTemperature(model, state)
        % Compute temperatures in each well
        
            % Get production well temperatures + injector/producer masks
            [Tprod, injector, producer] = model.getProductionWellTemperature(state);
            % Get injection well temperatures
            Tinj = model.getInjectionWellTemperature(state, Tprod);
            % Get well temperatures
            T = Tinj.*injector + Tprod.*producer;
            
        end
        
        %-----------------------------------------------------------------%
        function [T, injector, producer] = getProductionWellTemperature(model, state)
        % Get well temperature in production wells. This is computed based
        % on that the energy flux into all producing perforations of a well
        % equals the total energy production from that well. NOTE: This
        % funciton currently neglects that the density depends on pressure
        % and temperature
        
            % Get well mapping and perforation fluxes
            [map, q] = model.getProps(state, 'FacilityWellMapping', 'PhaseFlux');
            % Identify injecting and producing perforations
            qT = 0;
            for i = 1:numel(q)
                qT = qT + q{i};
            end
            T = zeros(numel(map.W), 1);
            injector = qT > 0;
            producer = ~injector;
            if ~any(producer)
                % No producing perforations
                injector = true(numel(map.W), 1);
                producer = ~injector;
                return;
            end
            % Get temperature in perforated cells
            Tperf = model.ReservoirModel.getProps(state, 'Temperature');
            Tperf = Tperf(map.cells);
            % Compute temperature in production wells
            qTprod = qT.*producer; % Only count fluxes into the well
            % Conservation of energy
            qMin    = -1e-16; % Avoid division by zero
            qTprodW = map.perforationSum*qTprod;
            Tprod   = (map.perforationSum*(Tperf.*qTprod))./min(qTprodW, qMin);
            % For wells producing less than qMin, we simply set the
            % production temperature equal to the mean temperature of the
            % preforated cells
            ok      = qTprodW < qMin;
            nc      = arrayfun(@(W) numel(W.cells), map.W);
            Tprod   = ok.*Tprod + ~ok.*map.perforationSum*Tperf./nc;
            % Identify injectors and producers
            q        = map.perforationSum*qT;
            injector = q > 0;
            producer = ~injector;
            % Tprod will have nans in wells with only injection perfs
            T           = model.AutoDiffBackend.convertToAD(zeros(size(value(Tprod))), Tprod);
            T(producer) = Tprod(producer);
            
        end
        
        %-----------------------------------------------------------------%
        function T = getInjectionWellTemperature(model, state, sample)
        % Get production well temperature
        
            % Get well mapping and perforation fluxes
            map = model.getProps(state, 'FacilityWellMapping');
            % Set injection temperature equal to W.T
            T = vertcat(map.W.T);
            % Make sure temperature is AD
            if ~isa(T, 'ADI')
                T = model.AutoDiffBackend.convertToAD(T, sample);
            end
            
        end
        
        function state = setWellValuesFromGroupControls(model, state0, state, dt, drivingForces)
        % Set well values for all wells that operates under group control
        
            % Check if groups are present
            groups = drivingForces.groups;
            ng     = numel(groups);
            if ng == 0, return; end
            % Compute well potential
            pot = model.computeWellPotential(state);
            % Loop through groups and update well controls
            for g = 1:ng
                group = groups{g};
                switch group.ctrlType
                    case 'rate'
                        state = model.setGroupWellRates(group, state, pot);
                    otherwise
                        error('Group control type not supported');
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [state, Tw] = setGroupWellRates(model, groupCtrl, state, pot)
        % Set well rates based on group target
        
            % Get well solutions
            ws = state.wellSol;
            map = model.getProps(state, 'FacilityWellMapping');
            if ~any(map.active), return; end
            % Get group wellSols
            mask = model.getGroupMask(state, groupCtrl.group);
            wsg  = ws(mask);
            % Get well rates and temperatures for all wells in group
            q = model.getWellRates(state);
            T = model.getWellTemperatures(state);
            if ischar(groupCtrl.ctrlVal)
                % Control value is the name of another group - this means
                % that we aim at a group rate equal to that groups rate,
                % with opposite sign
                maskp = model.getGroupMask(state, groupCtrl.ctrlVal);
                qtot  = -sum(q(maskp));
                Ttot  = sum(q(maskp).*T(maskp))./sum(q(maskp));
            elseif isnumeric(groupCtrl.ctrlVal)
                % Control value is a numeric variable
                qtot = groupCtrl.ctrlVal;
                Ttot = groupCtrl.T;
            end
            % Loop through wells in group
            active = true(numel(wsg),1);
            wix    = find(mask);
            Tw     = model.AutoDiffBackend.convertToAD(vertcat(ws.T), q);
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
        function pot = computeWellPotential(model, state)
        % Compute the well potential. This is the maximum
        % injection/production rate a well can manage when it is working at
        % its maximum/minimum bhp limit.
        
            map = model.getProps(state, 'FacilityWellMapping');
            for i = 1:numel(map.active)
                state.wellSol(map.active(i)).bhp = map.W(i).lims.bhp;
            end
            [pot, map] = model.getProps(state, 'PhaseFlux', 'FacilityWellMapping');
            pot = abs(map.perforationSum*sum(value(pot),2));
            
        end
        
        %-----------------------------------------------------------------%
        function mask = getGroupMask(model, state, name)
        % Get logical mask into wellSol for a given group
        
            map  = model.getProp(state, 'FacilityWellMapping');
            mask = arrayfun(@(W) strcmpi(W.group, name), map.W);
            
        end
        
        %-----------------------------------------------------------------%
        function q = getWellRates(model, state)
        % Get the total rate for each well
        
            nph = model.ReservoirModel.getNumberOfPhases();
            q = 0;
            % TODO: check which is faster - summing then slice, or slice
            % then sum
            for i = 1:nph
                q = q + state.FacilityState.surfacePhaseRates{i};
            end
            
        end
        
        %-----------------------------------------------------------------%
        function T = getWellTemperatures(model, state)
        % Get well temperature for each well
        
           isT = strcmpi(state.FacilityState.names, 'well_temperature');
           T   = state.FacilityState.primaryVariables{isT};
           
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Update state after convergence
            [state, report] = updateAfterConvergence@GenericFacilityModel(model, state0, state, dt, drivingForces);
            map = model.getProp(state, 'FacilityWellMapping');
            if strcmpi(model.thermalFormulation, 'none')
                Tw = num2cell(model.getWellTemperature(state));
                [state.wellSol(map.active).T] = deal(Tw{:});
            end
            if numel(map.active) < numel(state.wellSol)
                inactive = true(numel(state.wellSol),1);
                inactive(map.active) = false;
                if isfield(state0, 'wellSol')
                    [state.wellSol(inactive).T] = deal(state0.wellSol(inactive).T);
                end
            end
            
            % Compute produced/injected energy
            nw = numel(map.active);
            if nw > 0
                hAdv  = model.getProp(state, 'AdvectiveHeatFlux');
                hAdv  = expandMatrixToCell(hAdv);
                hCond = model.getProp(state, 'ConductiveHeatFlux');
                h     = hCond;
                for i = 1:model.ReservoirModel.getNumberOfPhases()
                    h = h + hAdv{i};
                end
                h = map.perforationSum*h;
                nwt = numel(state.wellSol);
                active = false(nwt, 1);
                active(map.active) = true;
                actIndex = zeros(nwt, 1);
                actIndex(map.active) = (1:numel(map.active))';
                for w = 1:nwt
                    if active(w)
                        hw = h(actIndex(w));
                    else
                        hw = 0;
                    end
                    state.wellSol(w).effect = hw;
                end
            end
            
        end
    end
    
end

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