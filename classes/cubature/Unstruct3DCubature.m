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
            cub.parentPos = (0:cub.numPoints:numParents*cub.numPoints)' + 1;
            
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
                  
                otherwise
                
                    error('Prescision not supported')
                
            end
            
        end
        
        function [x, w, n] = calculateWeights(cub, x)
            
            G = cub.G;
            
            n = size(x,1);
            
            if cub.prescision > 1
                
                dim = 3;
                cubTet = TetrahedronCubature(G, cub.prescision, cub.internalConn);
                basis = dgBasis(dim, cub.prescision, 'legendre');
                psi = basis.psi;
                nDof = basis.nDof;
                P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), nDof, nDof)';
                [W, xq, w, ii, jj, cellNo, faceNo] = cubTet.getCubature(1:G.cells.num, 'cell');

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