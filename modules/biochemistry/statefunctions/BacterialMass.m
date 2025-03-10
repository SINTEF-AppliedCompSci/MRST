classdef BacterialMass < StateFunction & ComponentProperty
    % Mass of each component, in each phase.
    properties (Access = protected)
        minimumDerivatives = 1.0e-3;
    end

    methods
        function gp = BacterialMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');
            gp.label = 'M_{bio}';
        end
        function mb = evaluateOnDomain(prop, model, state)

            s = model.getProps(state, 's');
            nbact = model.getProps(state, 'nbact');

            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            L_ix = model.getLiquidIndex();

            if iscell(s)
                Voln = s{L_ix};
            else
                Voln = s(:, L_ix);
            end
            Voln = max(Voln, 1.0e-8);
            mb = pv.*nbact.*Voln;
            mb = ensureMinimumDerivatives(prop,model, mb);

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

    function mass = ensureMinimumDerivatives(prop, model, mass)
        der = prop.minimumDerivatives;
        if isempty(der)
            return;
        end
        nc = numel(mass);
        if numel(der) == 1
            der = repmat(der, 1, nc);
        end
        isDiag = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
        if isDiag
            rowMajor = model.AutoDiffBackend.rowMajor;
            for c = 1:nc
                m = mass;
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
                mass = m;
            end
        else
            for c = 1:nc
                m = mass;
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
                mass = m;
            end
        end
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
