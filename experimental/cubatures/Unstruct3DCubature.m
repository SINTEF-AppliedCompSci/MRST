classdef Unstruct3DCubature < Cubature
    
    properties
        
        parentPos
        
    end
    
    methods
        
        function cub = Unstruct3DCubature(G, prescision, internalConn)
            
            cub = cub@Cubature(G, prescision, internalConn);
            [x, w, n] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            cub.dim = 3;
            
            numParents = G.cells.num;
            cub.pos = (0:cub.numPoints:numParents*cub.numPoints)' + 1;
            
        end
        
        function x = getCubaturePoints(cub)
            
            a = sqrt(3/5);
            b = sqrt(1/3);
            c = sqrt(3/7 + 2/7*sqrt(6/5));
            
            switch cub.prescision
            
                case 1
                
                    x = [0,0,0];
            
                case 2

                    x = [-a -b -a;
                          a -a -a;
                         -a  a -a;
                          a  a -a;
                         -a -a  a;
                          a -a  a;
                         -a  a  a;
                          a  a  a;
                          c  b  0;
                          0  0  0];
                  
%                 otherwise
%                 
%                     error('Prescision not supported')
                
            end
            
            if cub.prescision > 1
                G1 = computeGeometry(cartGrid([1,1,1], [2,2,2]));
                G1.nodes.coords = G1.nodes.coords - 1;
                G1 = computeVEMGeometry(G1);
                G1 = computeCellDimensions(G1); 
                cubTet = TetrahedronCubature(G1, cub.prescision, cub.internalConn);
                [~, x, ~, cellNo, ~] = cubTet.getCubature(1, 'volume');
                x = cubTet.transformCoords(x, cellNo);
                x = unique(x, 'rows');
                basis = dgBasis(3, cub.prescision, 'legendre');
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
            
            if cub.prescision > 1
                
                dim = 3;
                if isfield(G, 'parent')
                    knownCub = CoarseGrid3DCubature(G, cub.prescision, cub.internalConn);
                else
                    knownCub = TetrahedronCubature(G, cub.prescision, cub.internalConn);
                end
                
%                 knownCub = TetrahedronCubature(G, cub.prescision, cub.internalConn);
                basis = dgBasis(dim, cub.prescision, 'legendre');
                psi = basis.psi;
                nDof = basis.nDof;
                P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), nDof, nDof)';
                [~, xq, w, cellNo, ~] = knownCub.getCubature(1:G.cells.num, 'volume');

                xq = cub.transformCoords(xq, cellNo);
                
                
                I = cellfun(@(p) accumarray(cellNo, w.*p(xq)), psi, 'unif', false);

                rhs = zeros(nDof, G.cells.num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    i = I{dofNo};
                    i(abs(i) < tol) = 0;
                    rhs(dofNo, :) = i;
                end

                w = reshape(P\rhs, [], 1);
                
            else
                w = G.cells.volumes;
                
            end
            
            cellNo = reshape(repmat((1:G.cells.num), n, 1), [], 1);
            x = repmat(x, G.cells.num, 1);
            x = cub.transformCoords(x, cellNo, true);
            
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
