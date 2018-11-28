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
            % Construct cubature position vector
            cubature.pos = (0:cubature.numPoints:G.cells.num*cubature.numPoints)' + 1;
            
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
                cubTet = TetrahedronCubature(G, cubature.prescision, cubature.internalConn);
                [W, xq, wTet, cellNo, faceNo] = cubTet.getCubature((1:G.cells.num)', 'volume');
                xq     = cubature.transformCoords(xq, cellNo);
                % Moments
                M = cellfun(@(p) accumarray(cellNo, wTet.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, G.cells.num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                end
                % Matrix of basis functions evalauted at current quadrature
                % points
                P      = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), [], nDof)';
                % Compute significance
                significance = sum(P.^2,1);
                % Compute weights
                w      = reshape((P'/(P*P'))*rhs, [], 1);
                if cubature.reduce
                    % Try to eliminate least significant point until we have
                    % exactly nDof quadrature points
                    k = n;
                    tol = 1e-5;
                    if cubature.reduce
                        while k > nDof
                            [~, ix] = sort(significance);
                            xPrev = x;
                            wPrev = w;
                            for m = 1:numel(ix)
                                x(ix(m),:) = [];
                                P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), [], nDof)';
                                A = (P*P');
                                w = reshape((P'/A)*rhs, [], 1);
                                if all(isfinite(w)) && rcond(A) > tol
                                    significance = sum(P.^2,1);
                                    k = k-1;
                                    break
                                else
                                    x = xPrev;
                                    w = wPrev;
                                end
                            end

                        end
                    end
                    n = k;
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
            cellNo = reshape(repmat((1:G.cells.num), n, 1), [], 1);
            x = repmat(x, G.cells.num, 1);
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