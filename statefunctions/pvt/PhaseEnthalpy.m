classdef PhaseEnthalpy < StateFunction
%State function for computing the phase enthalpy

    methods
        %-----------------------------------------------------------------%
        function gp = PhaseEnthalpy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'enthalpy'}, 'state');
            if model.getNumberOfPhases == 2
                gp = gp.dependsOn({'pressure'}, 'state');
            end
            gp.label = 'h_{\alpha}';
        end
        
        %-----------------------------------------------------------------%
        function hph = evaluateOnDomain(prop, model, state)
            h   = model.getProps(state, 'enthalpy');
            nph = model.getNumberOfPhases();
            hph = cell(1, nph);
            [hph{:}] = deal(h);
            if nph == 2
                twoPhase = state.flag == 3;
                if any(twoPhase)
                    p      = model.getProp(state, 'pressure');
                    phases = model.getPhaseNames();
                    for i = 1:nph
                        hphi = prop.evaluateFluid(model, ['h', phases(i)], p);
                        hph{i}(twoPhase) = hphi(twoPhase);
                    end
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