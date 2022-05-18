classdef MolecularDiffusivity < StateFunction
%State function for computing molecular diffusivity
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function md = MolecularDiffusivity(model, varargin)
            md@StateFunction(model, varargin{:});
            md = md.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');
            md = md.dependsOn({'pressure', 'temperature'}, 'state');
            md.label = 'd_i';
        end
        
        %-----------------------------------------------------------------%
        function d = evaluateOnDomain(prop, model, state)
            % Get dependencies
            [p, T, pv] = model.getProps(state, 'pressure', 'temperature', 'PoreVolume');
            % Compute tourtuosity and porosity
            tau  = model.rock.tau;
            poro = pv./model.operators.vol;
            % Comoute diffusivity
            ncomp = model.getNumberOfComponents();
            d     = cell(ncomp,1);
            for i = 1:ncomp
                if model.Components{i}.molecularDiffusivity > 0
                    d{i} = model.Components{i}.molecularDiffusivity.*tau.*poro;
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