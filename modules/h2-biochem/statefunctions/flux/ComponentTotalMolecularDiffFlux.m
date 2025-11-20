classdef ComponentTotalMolecularDiffFlux < StateFunction
    % ComponentTotalMolecularDiffFlux - Total molecular diffusion flux
    % of components summed over all phases.

    methods
        function gp = ComponentTotalMolecularDiffFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('ComponentPhaseMolecularDiffFlux');
            gp.label = 'J_i^{mol,diff}';
        end

        function v = evaluateOnDomain(prop, model, state)
            % Get phase-wise diffusion fluxes
            J_phase = prop.getEvaluatedDependencies(state, 'ComponentPhaseMolecularDiffFlux');

            ncomp = model.getNumberOfComponents();
            v = cell(ncomp, 1);

            % Sum over phases for each component
            for c = 1:ncomp
                v{c} = 0;
                for ph = 1:model.getNumberOfPhases()
                    if ~isempty(J_phase{c, ph})
                        v{c} = v{c} + J_phase{c, ph};
                    end
                end
            end
        end
    end
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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