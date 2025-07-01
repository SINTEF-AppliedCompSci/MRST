classdef DiffusiveBactFlux < StateFunction & UpwindProperty
    % bacterial diffusive flux
    properties (Access = protected)
        Db = 10^(-8)*meter/second
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        function gp = DiffusiveBactFlux(model, varargin)                                  
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            gp@StateFunction(model);
            gp@UpwindProperty(upstr);
            gp.upwind_name = upwind_name;
            gp = gp.dependsOn(upwind_name);
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'Density'}, 'PVTPropertyFunctions');
           gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'Flux_{bio}';
        end
        function fluxbact = evaluateOnDomain(prop, model, state)
            % Get dependencies
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            nbact = model.getProps(state, 'nbact');
            n0 = model.nbactMax;
                       
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');
            s = model.getProps(state, 's');                
            L_ix = model.getLiquidIndex();

                
            if iscell(s)
                sL = s{L_ix};
                rhoL = rho{L_ix};

            else
                sL = s(:,L_ix);
                rhoL = rho(:,L_ix);
            end
            if iscell(sL)
                Voln = sL{1}.*rhoL{1};
            else
                Voln = sL.*rhoL;
            end
            Diffb = prop.Db;
            gradn = model.operators.Grad(nbact);
            %fluxbact = -0.001*prop.faceUpstream(model, state, flag{L_ix}, Diffb.*Voln).*gradn;

            fluxbact = -1.00e-12.*gradn;

                
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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