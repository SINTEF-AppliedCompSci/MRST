classdef ComponentTotalDiffusiveFlux < StateFunction
%State function for computing total diffusive flux of salt components due
%to concentration differences

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function df = ComponentTotalDiffusiveFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = df.dependsOn({'MolecularTransmissibility'});
            df = df.dependsOn({'ComponentPhaseMassFractions', 'Density'}, 'PVTPropertyFunctions');
            df.label = 'D_i';
        end
        
        %-----------------------------------------------------------------%
        function D = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.getEvaluatedDependencies(state, 'MolecularTransmissibility');
            ix       = ~cellfun(@isempty, T);
            T(ix)    = cellfun(@(T) T(model.operators.internalConn), T(ix), 'UniformOutput', false);
            [X, rho] = model.getProps(state, 'ComponentPhaseMassFractions', ...
                                             'Density'                    );
            ncomp = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();
            D     = cell(ncomp,1);
            [D{:}] = deal(0);
            for ph = 1:nph
                rhof = model.operators.faceAvg(rho{ph});
                for c = 1:ncomp
                    if ~any(value(T{c}))
                        continue
                    end
                    D{c} = D{c} - rhof.*T{c}.*op.Grad(X{c, ph});
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