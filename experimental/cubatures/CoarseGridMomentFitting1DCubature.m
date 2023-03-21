classdef CoarseGridMomentFitting1DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        
        reduce
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = CoarseGridMomentFitting1DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.reduce = false;
            cubature.dim       = 1;
            % Make cubature points and weights
            [x, w, n] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            % Construct cubature position vector
            cubature.pos = (0:cubature.numPoints:cubature.G.faces.num*cubature.numPoints)' + 1;
            
        end
           
        %-----------------------------------------------------------------%
        function [x, w, n] = makeCubature(cubature)
            
            % Dimension of cubature
            dim    = 1;
            G      = cubature.G;
            % Basis functions used in moment-fitting
            basis   = dgBasis(1, cubature.prescision, 'legendre');
            gradPsi = basis.grad_psi;
            psi     = basis.psi;
            nDof    = basis.nDof;

            % The starting point is a quadrature for a reference square
            % with more than nDof points
            k = cubature.prescision;
            nS = 0;
            [x, ~, nS, xR] = getLineCubaturePointsAndWeights(k);
            while nS < nDof
                k = k+1;
                [x, ~, nS, xR] = getLineCubaturePointsAndWeights(k);
            end
            
            
            cub = LineCubature(G, k, cubature.internalConn);
            % We use known cubature to calculate the moments
            parentCub = LineCubature(G.parent, cubature.prescision, cubature.internalConn);
%             [~, x, w, ~, faceNo] = cub.getCubature((1:G.faces.num)', 'face');
            
            faces = G.faces.fconn;
            
            [~, xp, wp, ~, fineFaceNo] = parentCub.getCubature(faces, 'face');
            % Map cubature points to reference coordinates
            
%             xp = (xp - G.parent.faces.centroids(fineFaceNo,:))./(G.parent.faces.areas(fineFaceNo)/2);
            
            
            
            ix = G.parent.faces.nodes(mcolon(G.parent.faces.nodePos(fineFaceNo), ...
                                             G.parent.faces.nodePos(fineFaceNo+1)-1));
            xf = G.parent.nodes.coords(ix,:);
            x0 = xf(1:2:end,:);
%             x1 = xf(2:2:end,:);
            
            xp = sqrt(sum((xp - x0).^2,2))./G.parent.faces.areas(fineFaceNo).*(xR(2) - xR(1)) + xR(1);
            
            faceNof   = rldecode((1:G.faces.num)', diff(G.faces.connPos), 1);
            np  = diff(parentCub.pos);            
            np  = accumarray(faceNof, np(faces));
            
%             cellNo = 
%             x     = cubature.transformCoords(x, faceNo);
%             xp    = cubature.transformCoords(xp, fineFaceNo);
            num   = G.faces.num;
            sgn = fineToCoarseSign(G);
            normals = G.parent.faces.normals./G.parent.faces.areas;
            normals = normals(fineFaceNo,:).*sgn;
            % Moments
            M = cellfun(@(p) accumarray(faceNof, wp.*p(xp)), psi, 'unif', false);
%             M{1} = G.faces.areas;
            % Compute right-hand side
            rhs = zeros(nDof, num);
            tol = eps(mean(G.cells.volumes));
            for dofNo = 1:nDof
                m = M{dofNo};
                m(abs(m) < tol) = 0;
                rhs(dofNo, :) = m;
            end
            moments = rhs(:);
            
            normals = G.faces.normals./G.faces.areas;
            [x,w,n] = fitMoments(x, basis, moments, num, 'reduce', false);
            
            x = cubature.mapCoords(x, xR);
%                 w = w;

%             % Map from reference to physical coordinates
%             if strcmp(type, 'face')
%                 % Face coordinates
%                 faceNo = reshape(repmat((1:G.faces.num), n, 1), [], 1);
%                 vec1   = G.faces.coordSys{1}(faceNo,:);
%                 vec2   = G.faces.coordSys{2}(faceNo,:);
%                 x = repmat(x, G.faces.num, 1).*(G.faces.dx(faceNo,:)/2);
%                 x = x(:,1).*vec1 + x(:,2).*vec2;
%                 x = x + G.faces.centroids(faceNo,:);
%             else
%                 % Cell coordinates
%                 cellNo = reshape(repmat((1:G.cells.num), n, 1), [], 1);
%                 x = repmat(x, G.cells.num, 1);
%                 x = cubature.transformCoords(x, cellNo, true);
%             end
            
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
            % Create mapping
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(cubature, x, faces, inverse)
            
            % Transform to/from reference coordinates
            G = cubature.G.parent;
            translation = -G.faces.centroids(faces,:);
            if isfield(G.faces, 'dx')
                scaling = 1./(G.faces.areas(faces)/2);
            else
                scaling = 1./(G.faces.diameters(faces)/(2*sqrt(G.griddim-1)));
            end
            
            if nargin < 4 || ~inverse
                xhat = (x + translation).*scaling;
            else
                xhat = x./scaling - translation;
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
