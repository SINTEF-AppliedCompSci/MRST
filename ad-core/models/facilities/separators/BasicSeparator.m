classdef BasicSeparator
    properties
        pressure
        T
    end
    
    methods
        function sep = BasicSeparator(varargin)
            sep = merge_options(sep, varargin{:});
        end
        
        function [phaseMassStreams, massFractions, densities] = separateComponentMassStream(sep, model, massStream)
            assert(iscell(massStream));
            ncomp = numel(massStream);
            nph = model.getNumberOfPhases();
            % Total mass stream of each phase
            phaseMassStreams = cell(1, nph);
            % Mass density of each phase
            densities = cell(1, nph);
            [phaseMassStreams{:}, densities{:}] = deal(0);
            % Mass fractions
            massFractions = cell(1, nph);
            [massFractions{:}] = deal(cell(1, ncomp));
            % Compute phase split based on simplified assumptions
            % (components act independently of each other)
            for c = 1:ncomp
                comp = model.Components{c};
                frac = comp.getPhaseCompositionSurface(model, [], sep.pressure, sep.T);
                for ph = 1:nph
                    result = massStream{c}*frac{ph};
                    massFractions{ph}{c} = result;
                    phaseMassStreams{ph} = phaseMassStreams{ph} + result;
                end
            end
            % We can now use a simple linear mixing rule to get densities
            for c = 1:ncomp
                rho = comp.getPurePhaseDensitySurface(model, [], sep.pressure, sep.T);
                for ph = 1:nph
                    densities{ph} = densities{ph} + massFractions{ph}{c}.*rho{ph};
                end
            end
        end
        
        function [moleStreams, moleFractions, molardensities] = separateComponentMoleStream(sep, model, molestream)
            assert(false, 'Not implemented in class');
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
