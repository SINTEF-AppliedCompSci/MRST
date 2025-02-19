classdef BactPermeability < StateFunction
%State function for computing pressure- and bacteria-dependent
%permeability

   properties
   end
   
   methods
       
        function perm = BactPermeability(model)
            perm@StateFunction(model);
            if model.dynamicFlowTrans
                perm = perm.dependsOn('pressure', 'state');
                if model.bacteriamodel                
                    perm = perm.dependsOn('nbact', 'state');
                end
            end
            perm.label = 'K';
        end
       
        function perm = evaluateOnDomain(prop, model, state)
            perm = model.rock.perm;
            if model.dynamicFlowTrans
                if model.bacteriamodel
                    [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                    perm = perm(p,nbact);
                else
                    p = model.getProps(state, 'pressure');
                    perm = perm(p,0);

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