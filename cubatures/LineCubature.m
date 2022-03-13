classdef LineCubature < Cubature
    % 1D line cubatureature for faces in 2D MRST grids
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = LineCubature(G, prescision)
            % Set up line cubatureature
            
            % Make sure we have a 2D grid structure
            assert(G.griddim == 2, 'LineCubature only supported for 2D grids')
            % Most of the construction is handled by parent class
            cubature     = cubature@Cubature(G, prescision, 1);
            cubature.dim = 1;
            
        end
        
        %-----------------------------------------------------------------%
        function x = mapCoords(cubature, x, xR)
            % Map cubatureature points from reference to physical coordinates
            
            G     = cubature.G;
            % Number of line vertices
            nPts  = 2;
            % Total number of cubatureature points
            nq    = size(x,1);
            nodes = G.faces.nodes(:,1);
            % Total number of lines
            nLin  = G.faces.num;
            % Node coordinates
            xn = G.nodes.coords(nodes,:);
            x1 = xn(1:2:end,:);
            x2 = xn(2:2:end,:);
            
            vec = rldecode(x2-x1, nq, 1);
            x1  = rldecode(x1, nq, 1);
            
            x = (repmat(x,nLin,1)-xR(1))./(xR(2) - xR(1)).*vec + x1;
            
        end
        
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature)
            % Make cubatureature
            
            G = cubature.G;
            % Total number of lines
            nLin = G.faces.num;
            % Get points and weights
            [x, w, n, xR] = getLineCubaturePointsAndWeights(cubature.prescision);
            % Map to physical coordinates
            x = cubature.mapCoords(x, xR);
            % Multiply weights by line lenghts
            linNo = reshape(repmat(1:G.faces.num, n, 1), [], 1);
            w = repmat(w, nLin, 1).*G.faces.areas(linNo);
            n = repmat(n, G.faces.num, 1);
            
            cubature = cubature.assignCubaturePointsAndWeights(x,w,n);
            
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
