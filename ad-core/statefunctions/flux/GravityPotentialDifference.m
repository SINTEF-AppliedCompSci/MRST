classdef GravityPotentialDifference < StateFunction
    % Difference in phase potential over a face due to gravity
    properties
        saturationWeighting = false; % Use saturation-weighted density for average
        weight = []; % Optional weighting matrix for gravity
    end
    
    methods
        function gp = GravityPotentialDifference(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            if gp.saturationWeighting
                gp = gp.dependsOn('s', 'state');
            end
            gp.label = 'g \rho_\alpha \Delta z';
        end
        function gRhoDz = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            nf = size(model.operators.N, 1);
            if norm(model.gravity) > 0 && nf > 0
                gdz = -model.getGravityGradient();
                nm = model.getPhaseNames();
                rho = prop.getEvaluatedExternals(model, state, 'Density');
                rho = expandMatrixToCell(rho);
                avg = model.operators.faceAvg;
                for i = 1:nph
                    if prop.saturationWeighting
                        s = model.getProp(state, ['s', nm(i)]);
                        rhof = avg(s.*rho{i})./max(avg(s), 1e-8);
                    else
                        rhof = avg(rho{i});
                    end
                    gRhoDz{i} = rhof.*gdz;
                end
            else
                [gRhoDz{:}] = deal(zeros(nf, 1));
            end
            w = prop.weight;
            if ~isempty(w)
                for i = 1:numel(gRhoDz)
                    gRhoDz{i} = w*gRhoDz{i};
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
