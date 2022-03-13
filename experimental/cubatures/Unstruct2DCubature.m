classdef Unstruct2DCubature < Cubature
    
    properties
        
        parentPos
        disc
        
    end
    
    methods
        
        function cub = Unstruct2DCubature(G, prescision, internalConn)
            
            cub = cub@Cubature(G, prescision, internalConn);
            cub.dim = 2;
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            
            if G.griddim == 2
                numParents = G.cells.num;
            else
                numParents = G.faces.num;
            end
            cub.pos = (0:cub.numPoints:numParents*cub.numPoints)' + 1;

        end
        
        function x = getCubaturePoints(cub)
            
            a = sqrt(3/5);
            b = sqrt(1/3);
            c = sqrt(3/7 - 2/7*sqrt(6/5));
            
            switch cub.prescision
                case 1
                    x = zeros(1, 2);
                case 2
                    x = [ 0  0;
                          0  b;
                         -a -a;
                          a -a;
                         -a  a;
                          a  a];
%                 otherwise
%                     error('Prescision not supported')                      
            end
                      
            if cub.prescision > 1
                
                G1 = computeGeometry(cartGrid([1,1], [2,2]));
                G1.nodes.coords = G1.nodes.coords - 1;
                G1 = computeVEMGeometry(G1);
                G1 = computeCellDimensions(G1);
                cubTri = TriangleCubature(G1, cub.prescision, cub.internalConn);
                [~, x, ~, cellNo, ~] = cubTri.getCubature(1, 'volume');
                x = cubTri.transformCoords(x, cellNo);
                x = unique(x, 'rows');
                basis = dgBasis(2, cub.prescision, 'legendre');
                nDof  = basis.nDof;
                psi = basis.psi;
                P = zeros(nDof, nDof);
                while rank(P) < nDof && cond(P) > 1
                    ix = randperm(size(x,1));
                    ix = ix(ix(1:nDof));
                    P  = reshape(cell2mat(cellfun(@(p) p(x(ix,:)), psi, 'unif', false)), nDof, nDof)';
                end
                x = x(ix,:);
            end
                    
                
        end
        
        function [x, w, n] = calculateWeights(cub, x)
            
            G = cub.G;
            
            n = size(x,1);
            
            if G.griddim > cub.dim
                type = 'face';
                elements = 1:G.faces.num;
            else
                type = 'volume';
                elements = 1:G.cells.num;
            end
            
            if cub.prescision > 1
                
                dim = 2;
                if 0
                if isfield(G, 'parent')
                    knownCub = CoarseGrid2DCubature(G, cub.prescision, cub.internalConn);
                else
                    knownCub = TriangleCubature(G, cub.prescision, cub.internalConn);
                end
                end
                knownCub = TriangleCubature(G, cub.prescision, cub.internalConn);
                basis = dgBasis(dim, cub.prescision, 'legendre');
                psi = basis.psi;
                nDof = basis.nDof;
                P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), nDof, nDof)';

                [W, xq, w, cellNo, faceNo] = knownCub.getCubature(elements, type);
                
                if G.griddim == 3
                    xq = cub.map2face(xq, faceNo);
                    count = faceNo;
                    num = G.faces.num;
                else
                    xq = cub.transformCoords(xq, cellNo);
                    count = cellNo;
                    num = G.cells.num;
                end
                
                I = cellfun(@(p) accumarray(count, w.*p(xq)), psi, 'unif', false);

                rhs = zeros(nDof, num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    i = I{dofNo};
                    i(abs(i) < tol) = 0;
                    rhs(dofNo, :) = i;
                end

                w = reshape(P\rhs, [], 1);
                
            else
                
                if G.griddim == 2
                    w = G.cells.volumes;
                else
                    w = G.faces.areas;
                end
                
            end
            
            if strcmp(type, 'face')
                
                faceNo = reshape(repmat((1:G.faces.num), n, 1), [], 1);
                vec1 = G.faces.coordSys{1}(faceNo,:);
                vec2 = G.faces.coordSys{2}(faceNo,:);
            
                x = repmat(x, G.faces.num, 1).*(G.faces.dx(faceNo,:)/2);
                x = x(:,1).*vec1 + x(:,2).*vec2;
                
                x = x + G.faces.centroids(faceNo,:);
                
            else
                
                cellNo = reshape(repmat((1:G.cells.num), n, 1), [], 1);
                x = repmat(x, G.cells.num, 1);
%                 xc     = G.cells.centroids(cellNo,:);
                x = cub.transformCoords(x, cellNo, true);
%                 dx = G.cells.dx(cellNo,:)/2;

%                 x = x.*dx + xc;
%                 
            end
            
%             x = cub.transformCoords(x, cellNo, true);
            
        end
        
        function x = map2face(cub, x, faceNo)
            
            G = cub.G;
            vec1 = G.faces.coordSys{1}(faceNo,:);
            vec2 = G.faces.coordSys{2}(faceNo,:);
            
            x = x - G.faces.centroids(faceNo,:);
            
%             nfe   = diff(G.faces.edgePos);
% %             nfe   = nfe(faceNo);
%             edges = G.faces.edges(cumsum(nfe));
%             edges = edges(faceNo);
%             nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
%     
%             vec1 = G.nodes.coords(nodes(2:2:end),:) - G.nodes.coords(nodes(1:2:end-1),:);
%             vec1 = vec1./sqrt(sum(vec1.^2,2));
%     
%             n    = G.faces.normals(faceNo,:)./G.faces.areas(faceNo);
%             vec2 = cross(vec1, n, 2);
            
            x = [sum(x.*vec1,2), sum(x.*vec2, 2)];
            x = x./(G.faces.dx(faceNo,:)/2);
            
%             x = (x - G.faces.centroids(faceNo,:))./(G.faces.dx(faceNo,:)/2);

            %             

%             V(:,1:2:2*numel(faceNo)) = reshape(vec1', 3, []);
%             V(:,2:2:2*numel(faceNo)) = reshape(vec2', 3, []);
%             
%             xx = x*V;
%             x = zeros(nq*numel(faceNo), G.griddim);
%             for dNo = 1:G.griddim
% 
%                 x(:, dNo) = reshape(xx(:, dNo:G.griddim:end), [], 1);
% 
%             end

%             faceNo = reshape(repmat(faceNo', nq, 1), [], 1);
%             cellNo = rldecode((1:G.cells.num)', (ncf - ncbf)*nq, 1);
%             nq = sum((1:G.cells.num)' == cellNo', 2);
%     
%             xf = disc.transformCoords(G.faces.centroids(faceNo,:), cellNo);
%             x = x + xf;
            
        end
        
        function [xhat, translation, scaling] = transformCoords(cub, x, cells, inverse)
            
            G = cub.G;
            translation = -G.cells.centroids(cells,:);
            if isfield(G.cells, 'dx')
                scaling = 1./(G.cells.dx(cells,:)/2);
            else
                scaling = 1./(G.cells.diameters(cells)/(2*sqrt(G.griddim)));
            end
            
            if nargin < 4 || ~inverse
                xhat = (x + translation).*scaling;
            else
                xhat = x./scaling - translation;
            end
               
        end
        
        function [x, w, n, cNo] = makeCubature(cub)
            
            x = cub.getCubaturePoints();
            [x, w, n] = cub.calculateWeights(x);
            
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
