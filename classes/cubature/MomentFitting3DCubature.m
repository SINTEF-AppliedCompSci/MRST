classdef MomentFitting3DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        
        reduce
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = MomentFitting3DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.reduce = false;
            % Make cubature points and weights
            [x, w, n] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            cubature.dim = 3;
            cubature.pos = [0; cumsum(n)] + 1;
            % Construct cubature position vector
%             cubature.pos = cumsum(
%             cubature.pos = (0:cubature.numPoints:G.cells.num*cubature.numPoints)' + 1;
            
        end
           
        %-----------------------------------------------------------------%
        function [x, w, n] = makeCubature(cubature)
            
            % Dimension of cubature
            dim    = 3;
            % Basis functions used in moment-fitting
            basis  = dgBasis(dim, cubature.prescision, 'legendre');
            psi    = basis.psi;
            nDof   = basis.nDof;
            % The starting point is a quadrature for a reference square
            % with more than nDof points
            k = 0;
            nS = 0;
            while nS < nDof
                [x, ~, nS] = getCubeCubaturePointsAndWeights(cubature.prescision + k);
                k = k+1;
            end
            n = size(x,1);
            G = cubature.G;
            
            if cubature.prescision > 1
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % We use known cubature to calculate the moments
                                % We use known cubature to calculate the moments
                if isfield(G, 'parent')
                    knownCub = CoarseGrid3DCubature(G, cubature.prescision, cubature.internalConn);
                else
                    knownCub = TetrahedronCubature(G, cubature.prescision, cubature.internalConn);
                end
%                 knownCub = TetrahedronCubature(G, cubature.prescision, cubature.internalConn);
                [~, xq, wq, cellNo] = knownCub.getCubature((1:G.cells.num)', 'volume');
                xq     = cubature.transformCoords(xq, cellNo);
                % Moments
                M = cellfun(@(p) accumarray(cellNo, wq.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, G.cells.num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                    M{dofNo} = m;
                end
                [x,w,n] = fitMoments2(x, basis, M, 'equal', G.equal);
                
                if numel(w) == 1
                    w = repmat(w{:}, G.cells.num, 1);
                    x = repmat(x{:}, G.cells.num, 1);
                    n = n*ones(G.cells.num,1);
                else
                    w = vertcat(w{:});
                    x = vertcat(x{:});
                end
            else
                n = 1;
                % Precision is 1, use midpoint rule
                if G.griddim == 2
                    w = G.cells.volumes;
                else
                    w = G.faces.areas;
                end
            end
            
            % Map from reference to physical coordinates
            cellNo = rldecode((1:G.cells.num)', n, 1);
            x = cubature.transformCoords(x, cellNo, true);
            
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(cubature, x, cells, inverse)
            % Transform to/from reference coordinates
            G = cubature.G;
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
        
    end
    
end