classdef AquiferModel < PhysicalModel
    
    properties
        aquifers % aquifers data table
        aquind % structure which explains the fields in the aquifers data table
        aquiferprops
        initvals
        ReservoirModel % reservoir model the aquifer belongs to
    end
    
    methods
        function model = AquiferModel(reservoirModel, aquifers, aquind, ...
                                      aquiferprops, initvals, varargin)
            model = model@PhysicalModel([]);
            model.ReservoirModel = reservoirModel;
            model.aquifers     = aquifers;
            model.aquind       = aquind;
            model.aquiferprops = aquiferprops;
            model.initvals     = initvals;
        end
        
        function state = validateState(model, state)
        if ~isfield(state, 'aquiferSol')
            state = initStateAquifer(model, state);
        end
        end
        
        % -----------------------------------------------------------------
        
        function model = validateModel(model, reservoirModel)
            model.ReservoirModel = reservoirModel;
        end
        
        % -----------------------------------------------------------------
        
        function state = initStateAquifer(model, state)
        p   = model.initvals.pressures;
        vol = model.initvals.volumes;
        
        naq = numel(p);
        sol = repmat(struct('pressure', [], 'volume', [], 'flux', []), [naq, 1]);
        for k = 1:naq
            [sol(k).pressure, sol(k).volume] = deal(p(k), vol(k));
        end
        state.aquiferSol = sol;
        end
        
        % -----------------------------------------------------------------
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            report = [];
            p0   = vertcat(state.aquiferSol.pressure);
            vol0 = vertcat(state.aquiferSol.volume);
            
            q = computeAquiferFluxes(model, state, dt);
            
            aq    = model.aquifers;
            ix    = model.aquind;
            props = model.aquiferprops;
            
            conn = aq(:, ix.conn);
            aquid = aq(:, ix.aquid);
            nconn = size(conn, 1);
            naq = max(aquid);
            aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
            
            Q = aquid2conn'*q;
            Q = dt.*Q;
            C = props.C;
            p = p0 - Q./(C.*vol0);
            vol = vol0 - Q;
            
            for k = 1:numel(p)
                state.aquiferSol(k).pressure = p(k);
                state.aquiferSol(k).volume   = vol(k);
                state.aquiferSol(k).flux     = q(logical(aquid2conn(:,k)));
            end
        end
        
        % -----------------------------------------------------------------
        
        function q = computeAquiferFluxes(model, state, dt)
            % compute aquifer connection fluxes
            rmodel = model.ReservoirModel;
                      
            p_aq = vertcat(state.aquiferSol.pressure);
            V_aq = vertcat(state.aquiferSol.volume);

            aq = model.aquifers;
            ix   = model.aquind;
            
            alpha     = aq(:, ix.alpha);
            J         = aq(:, ix.J);
            conn      = aq(:, ix.conn);
            depthconn = aq(:, ix.depthconn);
            depthaq   = aq(:, ix.depthaq);
            C         = aq(:, ix.C);
            aquid     = aq(:, ix.aquid);
            
            nconn = size(conn, 1);
            naq = max(aquid);
            aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
            
            p_aq = aquid2conn*p_aq;
            V_aq = aquid2conn*V_aq;
            
            p  = rmodel.getProps(state, 'PhasePressures');
            ph = rmodel.getPhaseNames();
            isw = ph == 'W';
            pW  = p{isw}(conn);

            b = rmodel.getProps(state, 'ShrinkageFactors');
            bW = b{isw}(conn);
            rhoW = bW.*rmodel.fluid.rhoWS;
            
            Tc = C.*V_aq./J;
            if dt == 0
                coef = 1;
            else
                coef = (1 - exp(-dt./Tc))./(dt./Tc);
            end
            
            g = rmodel.gravity(3);
            q = alpha.*J.*(p_aq - pW + rhoW.*g.*(depthconn - depthaq)).*coef;
            
        end

        % -----------------------------------------------------------------
        
        function eqs = addAquifersContribution(model, eqs, names, state,  dt)
            q = computeAquiferFluxes(model, state, dt);
            wind = strcmpi('water', names);
            scale = 1;
            if isa(model.ReservoirModel, 'GenericReservoirModel')
                rhoS  = model.ReservoirModel.getSurfaceDensities();
                scale = rhoS(wind);
            end
            ix   = model.aquind;
            conn = model.aquifers(:, ix.conn);
            eqs{wind}(conn) = eqs{wind}(conn) - q*scale;
        end
        
        % -----------------------------------------------------------------
        
        function state = updateAquiferFluxes(model, state, dt)
            if nargin < 3
                dt = 0;
            end
            q = computeAquiferFluxes(model, state, dt);
            [aq, ix] = deal(model.aquifers, model.aquind);
            aquid    = aq(:, ix.aquid);
            for k = 1:naq
                state.aquiferSol(k).flux = q(aquid==k);
            end
        end
    end
end

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
