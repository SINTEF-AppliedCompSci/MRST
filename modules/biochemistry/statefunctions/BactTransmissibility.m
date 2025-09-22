classdef BactTransmissibility < StateFunction
    %State function for computing bacterial transmissibility in the fluid or rock

    properties
        postfix = '';
    end

    methods
        %-----------------------------------------------------------------%
        function btrans = BactTransmissibility(model, postfix)
            btrans@StateFunction(model);
            btrans.postfix = postfix;
            if isfield(model.fluid, ['transMult', postfix])
                btrans = btrans.dependsOn({'pressure', 'T'}, 'state');
            end
            name        = strcat('T', postfix);
            btrans.label = ['T_{', postfix, '}'];
            assert(isfield(model.operators, name));
            T = value(model.operators.(name));
            assert(all(isfinite(T)))
            if any(T < 0)
                warning('Negative transmissibility in %d interfaces', sum(T < 0));
            end
            btrans.outputRange = [0, inf];
        end

        %-----------------------------------------------------------------%
        function Tb = evaluateOnDomain(trans, model, state)
            name  = strcat('T'        , trans.postfix);
            mname = strcat('transMult', trans.postfix);
            Tb = model.operators.(name);
            if isfield(model.fluid, mname)
                % Get pressure and temperature
                [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                p      = model.operators.faceAvg(p);
                nbact      = model.operators.faceAvg(nbact);
                % Apply multiplier
                transMult = trans.evaluateFluid(model, mname, p, nbact);
                Tb        = transMult.*Tb;
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