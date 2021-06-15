classdef EOSSeparator < BasicSeparator
    properties
        EOSModel
    end
    
    methods
        function sep = EOSSeparator(varargin)
            sep = merge_options(sep, varargin{:});
        end
        
        function [massstreams, massfractions, densities] = separateComponentMassStream(sep, model, molestream)
            assert(false, 'Use "moles" as mode for EOSSeparator');
        end
        
        function [molestreams, molefractions, densities] = separateComponentMoleStream(sep, model, molestream)
            if isempty(sep.EOSModel)
                eos = model.EOSModel;
            else
                eos = sep.EOSModel;
            end
            assert(iscell(molestream));
            totMole = 0;
            ncomp = numel(molestream);
            for i = 1:ncomp
                totMole = totMole + molestream{i};
            end
            z = molestream;
            for i = 1:ncomp
                z{i} = z{i}./totMole;
            end
            zD = ensureMinimumFraction(value(z), eos.minimumComposition);
            n = numelValue(z{1});
            p = repmat(sep.pressure, n, 1);
            temp = repmat(sep.T, n, 1);
            % Build state
            state = struct('pressure', p, 'T', temp, 'components', zD);
            % Perform flash
            solver = getDefaultFlashNonLinearSolver();
            state = eos.validateState(state);
            state = solver.solveTimestep(state, 1, eos);
            % Set derivatives, if we are using AD
            [Z_L, Z_V] = deal(state.Z_L, state.Z_V);
            useAD = isa(z{1}, 'ADI') || isa(p, 'ADI') || isa(temp, 'ADI');
            if useAD
                [x, y, L] = eos.getPhaseFractionAsADI(state, p, temp, z);
                acf = eos.CompositionalMixture.acentricFactors;
                [A_ij, Bi] = eos.getMixingParameters(p, temp, acf, true);
                [~, A_L, B_L] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
                [~, A_V, B_V] = eos.getPhaseMixCoefficients(y, A_ij, Bi);
                Z_L = eos.setZDerivatives(Z_L, A_L, B_L);
                Z_V = eos.setZDerivatives(Z_V, A_V, B_V);
            else
                x = state.x;
                y = state.y;
                L = state.L;
            end
            rhoL = eos.PropertyModel.computeMolarDensity(eos, p, x, Z_L, temp, true);
            rhoV = eos.PropertyModel.computeMolarDensity(eos, p, y, Z_V, temp, false);

            densities = {rhoL, rhoV};
            moleL = L.*totMole;
            moleV = (1-L).*totMole;
            
            molestreams = {moleL, moleV};
            if ~useAD
                x = expandMatrixToCell(x);
                y = expandMatrixToCell(y);
            end
            molefractions = {x, y};
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
