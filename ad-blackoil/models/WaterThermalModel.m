classdef WaterThermalModel < WaterModel
    % Single phase water model with thermal effects. Should be considered
    % experimental and intentionally undocumented, as this feature is
    % subject to change in the future.
    properties
        thermal = true;
    end
    
    methods
        function model = WaterThermalModel(G, rock, fluid, varargin)
            model = model@WaterModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
            
            rock_heat = struct('perm',rock.lambdaR);
            T_r=computeTrans(G,rock_heat);
            cf = G.cells.faces(:,1);
            nf = G.faces.num;
            T_r  = 1 ./ accumarray(cf, 1./T_r, [nf, 1]);
            model.operators.T_r_all=T_r;
            intInx = all(G.faces.neighbors ~= 0, 2);
            model.operators.T_r = T_r(intInx);
            % Setup operators
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ReservoirModel(model);
            %
            forces.bcT = [];
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsWaterThermal(state0, state, model,...
                dt, ...
                drivingForces,...
                varargin{:});
            
        end
        
        function rhoS = getSurfaceDensities(model)
            active = model.getActivePhases();
            props = {'rhoWS', 'rhoOS', 'rhoGS'};
            rhoS = cellfun(@(x) model.fluid.(x), props(active));
        end
        
        function [eqs, names, types, wellSol] = insertWellEquations(model, eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, hW,components, dt, opt)
            % Utility function for setting up the well equations and adding
            % source terms for black-oil like models. Note that this currently
            % assumes that the first nPh equations are the conservation
            % equations, according to the canonical MRST W-O-G ordering,
            fm = model.FacilityModel;
            active=model.getActivePhases;
            nPh = nnz(active);
            assert(numel(eqs) == nPh+1);
            dissolved =[];
            [eqs, names, types, wellSol, src] = insertWellEquations@ReservoirModel(model, eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, ...
                                                     wellMap, p, mob, rho, ...
                                                     dissolved, components, ...
                                                     dt, opt);
            rhoS = model.getSurfaceDensities();
            rhoS=rhoS(active);
            wc = fm.getWellCells();
            actw = fm.getWellStatusMask(wellSol);
            W = fm.getWellStruct(actw);
           
            % get the entalphy of the wells
            % should probably be moved into a separate well facility model
            nw = numel(W);
            nwc = nan(nw,1);
            hWW = nan(nw,1);
            for i = 1:nw
                nwc(i) = numel(W(i).cells);
                hWW(i) = W(i).hW;
            end
            hWW=rldecode(vertcat(hWW),nwc);
            % hard code well part for energy equation use upwind
            % assume no heat conduction
            % no contronling of heat conduction part
            % should be added if adjoint is needed
            cqs= src.phaseMass{1}./rhoS(1);
            hFw=rhoS(1)*hW(wc);
            hFww=rhoS(1).*hWW;%Rw*vertcat(W.hW);
            hFw(cqs>0)=hFww(cqs>0);
            eqs{2}(wc)= eqs{2}(wc) - hFw.*cqs;
            
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