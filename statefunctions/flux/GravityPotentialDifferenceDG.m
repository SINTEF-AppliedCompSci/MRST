classdef GravityPotentialDifferenceDG < GravityPotentialDifference
    
    methods
        function gp = GravityPotentialDifferenceDG(varargin)
            gp = gp@GravityPotentialDifference(varargin{:});
        end
        function gRhoDz = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            
            g = model.getGravityVector();
            x = model.G.cells.centroids(state.cells,:);
            gdz = model.operators.Grad(x)*g';
                
            if norm(model.gravity) > 0
                nm = model.getPhaseNames();
                rho = model.getProp(state, 'Density');
                avg = model.operators.faceAvg;
                for i = 1:nph
                    if prop.saturationWeighting
                        s = model.getProp(state, ['s', nm(i)]);
                        rhof = avg(s.*rho{i})./max(avg(s), 1e-8);
                    else
                        rhof = avg(rho{i});
                    end
                    gRhoDz{i} = - rhof.*gdz;
                end
            else
                [gRhoDz{:}] = deal(gdz);
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
