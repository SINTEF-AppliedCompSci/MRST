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
            dim = 3;
            G   = cubature.G;
            % Basis functions used in moment-fitting
            basis = dgBasis(dim, cubature.prescision, 'legendre');
            psi   = basis.psi;
            nDof  = basis.nDof;

            if cubature.prescision <=1
                % Precision is 1, use midpoint rule
                [x, ~, n] = getCubeCubaturePointsAndWeights(cubature.prescision);
                w = G.cells.volumes;
                x = repmat(x, G.cells.num,1);
                n = repmat(n, G.cells.num,1);
                cellNo = (1:G.cells.num)';
                
            else
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % The starting point is a quadrature for a reference square
                % with more than nDof points
                if 1
                    G1 = computeGeometry(cartGrid([1,1,1], [2,2,2]));
                    G1.nodes.coords = G1.nodes.coords - 1;
                    G1 = computeVEMGeometry(G1);
                    G1 = computeCellDimensions2(G1);
                    cubTet = TetrahedronCubature(G1, cubature.prescision, cubature.internalConn);
                    [~, x, ~, cellNo, ~] = cubTet.getCubature(1, 'volume');
                    x = cubTet.transformCoords(x, cellNo);
                    x = unique(x, 'rows');
                 else
                    k = 0;
                    nS = 0;
                    while nS < 2*nDof
                        [x, ~, nS] = getCubeCubaturePointsAndWeights(cubature.prescision + k);
                        k = k+1;
                    end
                end
                
                % We use known cubature to calculate the moments
                if isfield(G, 'parent')
                    knownCub = CoarseGrid3DCubature(G, cubature.prescision, cubature.internalConn);
                else
                    knownCub = TetrahedronCubature(G, cubature.prescision, cubature.internalConn);
                end
                [~, xq, wq, cellNo] = knownCub.getCubature((1:G.cells.num)', 'volume');
                wq = wq./G.cells.volumes(cellNo);
                xq     = cubature.transformCoords(xq, cellNo);
                % Moments
                M = cellfun(@(p) accumarray(cellNo, wq.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, G.cells.num);
                tol = eps(1);
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                    M{dofNo} = m;
                end
                [x,w,n] = fitMoments3(x, basis, M);
                cellNo = rldecode((1:G.cells.num)', n, 1);
                w = w.*G.cells.volumes(cellNo);
            end
            
            % Map from reference to physical coordinates
            x = cubature.transformCoords(x, cellNo, true);
            
        end        
    end
    
end