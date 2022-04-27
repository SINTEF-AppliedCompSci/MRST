classdef WellAdvectiveHeatFlux < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = WellAdvectiveHeatFlux(model, varargin)
            
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseFlux', 'FacilityWellMapping'});
            if model.implicitConnectionDP
                gp = gp.dependsOn('ConnectionPressureDrop');
            end
            gp.label = 'q_{avd}';
            
        end
        
        %-----------------------------------------------------------------%
        function q = evaluateOnDomain(prop, model, state)
            
            % Get well phase fluxes and well mapping
            [v, map] = prop.getEvaluatedDependencies(state, 'PhaseFlux', 'FacilityWellMapping');
            if model.implicitConnectionDP
                cdp = prop.getEvaluatedDependencies(state, 'ConnectionPressureDrop');
            else
                cdp = vertcat(state.wellSol(map.active).cdp);
            end
            % Get enthalpy and density in perforated cells
            [h, rho] = model.ReservoirModel.getProps(state, 'PhaseEnthalpy', 'Density');
            h   = cellfun(@(h)   h(map.cells)  , h  , 'UniformOutput', false);
            rho = cellfun(@(rho) rho(map.cells), rho, 'UniformOutput', false);
            % Identify injecting perforations
            vt = 0;
            for i = 1:numel(v)
                vt = vt + v{i};
            end
            injector = value(vt) > 0;
            if any(injector)
                % For injecting perforations, compute enthalpy as seen from
                % wellore and update advecitve heat flux accordingly
                % Get bottom-hole pressure from well variables
                isbhp = strcmpi(state.FacilityState.names, 'bhp');
                if any(isbhp)
                    bhp = state.FacilityState.primaryVariables{isbhp}(map.perf2well);
                else
                    bhp = vertcat(state.wellSol(map.active).bhp);
                    bhp = bhp(map.perf2well);
                end
                cp  = bhp + cdp;
                % Get temperature from well variables
                isT = strcmpi(state.FacilityState.names, 'well_temperature');
                if any(isT)
                    Tw = state.FacilityState.primaryVariables{isT}(map.perf2well);
                else
                    Tw = vertcat(map.W.T);
                    Tw = Tw(map.perf2well);
                end
                % Get NaCl concentraion if present
                xw     = vertcat(map.W.components);
                Xw     = model.ReservoirModel.getMassFraction(xw);
                Xw     = Xw(map.perf2well,:);
                cnames = model.ReservoirModel.getComponentNames();
                isNaCl = strcmpi(cnames, 'NaCl');
                if any(isNaCl), Xw = Xw(:,isNaCl); else, Xw = []; end
                % Compute enthalpy and density in wellbore
                phases = model.ReservoirModel.getPhaseNames();
                for i = 1:numel(phases)
                    ix = model.ReservoirModel.getPhaseIndex(phases(i));
                    hWell = feval(model.ReservoirModel.fluid.(['h', phases(i)]), cp, Tw);
                    h{ix}(injector) = hWell(injector);
                    rhoWell = feval(model.ReservoirModel.fluid.(['rho', phases(i)]), cp, Tw, Xw);
                    rho{ix}(injector) = rhoWell(injector);
                end
            end
            % Compute advective heat flux in/out of wellbore
            q = cellfun(@(rho,h,v) rho.*v.*h, rho, h, v, 'UniformOutput', false);
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