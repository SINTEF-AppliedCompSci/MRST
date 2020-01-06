classdef AquiferModel < PhysicalModel
    
    properties
        aquifers % aquifers data table
        aquind % structure which explains the fields in the aquifers data table
        aquiferprops
        reservoirModel % reservoir model the aquifer belongs to
    end
    
    methods
        function model = AquiferModel(reservoirModel, aquifers, aquind, ...
                                      aquiferprops, varargin)
            model = model@PhysicalModel([]);
            model.reservoirModel = reservoirModel;
            model.aquifers     = aquifers;
            model.aquind       = aquind;
            model.aquiferprops = aquiferprops;
        end
        
        function state = validateState(model, state)
        % Check parent class
            fn = model.getVariableField('aquiferfluxes');
            if ~isfield(state, fn)
                dt = 0;
                q = computeAquiferFluxes(model, state, dt);
                state.aquiferfluxes = q;
            end
            clear fn
        end
        
        % % --------------------------------------------------------------------%
        % function [fn, index] = getVariableField(model, name, varargin)
        %     switch(lower(name))
        %       case {'aquiferpressures', 'aquifervolumes'}
        %         fn = lower(name);
        %         index = 1;
        %       case 'aquiferfluxes'
        %         fn = lower(name);
        %         index = 1;
        %       otherwise
        %         % Basic phases are known to the base class
        %         [fn, index] = getVariableField@PhysicalModel(model, name, ...
        %                                                      varargin{:});
        %     end
        % end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            reservoirModel = model.reservoirModel;
            
            p_aq = model.getProp(state, 'aquiferpressures');
            V_aq = model.getProp(state, 'aquifervolumes');
            
            q = computeAquiferFluxes(model, state, dt);
            
            aquifers     = model.aquifers;
            aquind       = model.aquind;
            aquiferprops = model.aquiferprops;
            conn = aquifers(:, aquind.conn);
            aquid = aquifers(:, aquind.aquid);
            nconn = size(conn, 1);
            naq = max(aquid);
            aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
            
            Q = aquid2conn'*q;
            Q = dt.*Q;
            C = aquiferprops.C;
            p_aq = p_aq - Q./(C.*V_aq);
            V_aq = V_aq - Q;
            
            % TODO: reorder by aquifer
            state.aquifersol = struct('aquiferpressures', p_aq, ...
                                      'aquifervolumes'  , V_aq, ...
                                      'aquiferfluxes'   , q);
        end
        
        function q = computeAquiferFluxes(model, state, dt)

            reservoirModel = model.reservoirModel;
            [p, sW] = reservoirmodel.getProps(state, 'pressure', 'water');
            
            % TODO : fetch in right-way (fix below)
            p_aq = state.aquifersol.aquiferpressures;
            V_aq = state.aquifersol.aquifervolumes;

            aquifers = model.aquifers;
            aquind   = model.aquind;
            
            alpha     = aquifers(:, aquind.alpha);
            J         = aquifers(:, aquind.J);
            conn      = aquifers(:, aquind.conn);
            depthconn = aquifers(:, aquind.depthconn);
            depthaq   = aquifers(:, aquind.depthaq);
            C         = aquifers(:, aquind.C);
            aquid     = aquifers(:, aquind.aquid);
            
            nconn = size(conn, 1);
            naq = max(aquid);
            aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
            
            p_aq = aquid2conn*p_aq;
            V_aq = aquid2conn*V_aq;

            p  = p(conn);
            sW = sW(conn);
            
            fluid = model.fluid;
            pcOW = 0;
            ph = model.getPhaseNames();
            isw = (ph == 'W');
            if isfield(fluid, 'pcOW') && ~isempty(sW)
                pcOW = model.getProps(state, 'CapillaryPressure');
                pcOW = pcOW{isw}(conn);
            end

            b = model.getProps(state, 'ShrinkageFactors');
            bW = b{isw}(conn);
            rhoW = bW.*fluid.rhoWS;
            
            Tc = C.*V_aq./J;
            if dt == 0
                coef = 1;
            else
                coef = (1 - exp(-dt./Tc))./(dt./Tc);
            end
            
            g = model.gravity(3);
            q = alpha.*J.*(p_aq + pcOW - p + rhoW.*g.*(depthconn - depthaq)).*coef;
            
        end

        function eqs = addAquifersContribution(model, eqs, names, state,  dt)
            
            reservoirModel = model.reservoirModel;
           
            q = computeAquiferFluxes(model, state, dt);
            
            wind = strcmpi('water', names);
            aquifers = model.aquifers;
            aquind   = model.aquind;
            conn = aquifers(:, aquind.conn);
            eqs{wind}(conn) = eqs{wind}(conn) - q;

        end
        
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
