classdef BactPorosity < StateFunction
%State function for computing pressure- and bacteria-dependent
%permeability

   properties
   end
   
   methods
       
        function poro = BactPorosity(model)
            poro@StateFunction(model);
            if model.dynamicFlowTrans
                poro = poro.dependsOn({'pressure', 'nbact'}, 'state');
            end
            poro.label = 'Phi';
        end
       
       function poro = evaluateOnDomain(prop, model, state)
           poro = model.rock.poro;
           if model.dynamicFlowPv
               [p, nbact] = model.getProps(state, 'pressure', 'nbact');
               poro = poro(p,nbact);
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