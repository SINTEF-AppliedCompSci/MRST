classdef ComponentTotalMass <  StateFunction & ComponentProperty
    % The total mass of each component, given per cell
    properties (Access = protected)
        minimumDerivatives = [];
    end
    
    methods
        function gp = ComponentTotalMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('ComponentPhaseMass');
            gp.label = 'M_i';
        end

        function mass = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            mass = cell(ncomp, 1);
            phase_mass = prop.getEvaluatedDependencies(state, 'ComponentPhaseMass');
            for c = 1:ncomp
                % Loop over phases where the component may be present
                for ph = 1:nph
                    % Check if present
                    m = phase_mass{c, ph};
                    if ~isempty(m)
                        if isempty(mass{c})
                            mass{c} = m;
                        else
                            mass{c} = mass{c} + m;
                        end
                    end
                end
            end
            mass = prop.ensureMinimumDerivatives(mass);
        end
        
        function prop = setMinimumDerivatives(prop, der)
            if nargin < 2
                % Just set to a very small value, since nothing was
                % specified.
                der = 1e-10;
            end
            prop.minimumDerivatives = der;
        end
        
        function der = getMinimumDerivatives(prop)
            der = prop.minimumDerivatives;
        end

        function mass = ensureMinimumDerivatives(prop, mass)
            der = prop.minimumDerivatives;
            if isempty(der)
                return;
            end
            nc = numel(mass);
            if numel(der) == 1
                der = repmat(der, 1, nc);
            end
            isDiag = isa(prop.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiag
                rowMajor = prop.AutoDiffBackend.rowMajor;
                for c = 1:nc
                    m = mass{c};
                    if isnumeric(m) || size(m.jac{1}.diagonal, 2 - rowMajor) < c
                        continue
                    end
                    d = der(c);
                    if rowMajor
                        bad = abs(m.jac{1}.diagonal(c, :)) < d;
                        m.jac{1}.diagonal(c, bad) = d;
                    else
                        bad = abs(m.jac{1}.diagonal(:, c)) < d;
                        m.jac{1}.diagonal(bad, c) = d;
                    end
                    mass{c} = m;
                end
            else
                for c = 1:nc
                    m = mass{c};
                    d = der(c);
                    if isnumeric(m)
                        continue;
                    end
                    if numel(m.jac) < c
                        % Not matching number of derivatives - stop early
                        break;
                    end
                    J = m.jac{c};
                    [n, l] = size(J);
                    if n ~= l
                        continue;
                    end
                    diagonal = diag(J);
                    bad = abs(diagonal) < d;
                    if any(bad)
                        m.jac{c} = m.jac{c} + sparse(1:n, 1:n, bad.*d, n, n);
                    end
                    mass{c} = m;
                end
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
