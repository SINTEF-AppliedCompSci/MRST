classdef HeatTransmissibility < StateFunction
%State function for computing heat transmissibility in the fluid or rock
    
    properties
        postfix = '';
    end
    
    methods
        %-----------------------------------------------------------------%
        function trans = HeatTransmissibility(model, postfix)
            trans@StateFunction(model);
            trans.postfix = postfix;
            if isfield(model.fluid, ['transMult', postfix])
                trans = trans.dependsOn({'pressure', 'T'}, 'state');
            end
            name        = strcat('T', postfix);
            trans.label = ['T_{', postfix, '}'];
            assert(isfield(model.operators, name));
            T = value(model.operators.(name));
            assert(all(isfinite(T)))
            if any(T < 0)
                warning('Negative transmissibility in %d interfaces', sum(T < 0));
            end
            trans.outputRange = [0, inf];
        end

        %-----------------------------------------------------------------%
        function Th = evaluateOnDomain(trans, model, state)
            name  = strcat('T'        , trans.postfix);
            mname = strcat('transMult', trans.postfix);
            Th = model.operators.(name);
            if isfield(model.fluid, mname)
                % Get pressure and temperature
                [p, T] = model.getProps(state, 'pressure', 'T');
                p      = model.operators.faceAvg(p);
                T      = model.operators.faceAvg(T);
                % Apply multiplier
                transMult = trans.evaluateFluid(model, mname, p, T);
                Th        = transMult.*Th;
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