classdef MomentFitting2DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        
        reduce
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = MomentFitting2DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.reduce = false;
            cubature.dim       = 2;
            % Make cubature points and weights
            [x, w, n] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            % Construct cubature position vector
            if G.griddim == 2
                numParents = G.cells.num;
            else
                numParents = G.faces.num;
            end
            cubature.pos = [0; cumsum(n)] + 1;
%             cubature.pos = (0:cubature.numPoints:numParents*cubature.numPoints)' + 1;
            
        end
           
        %-----------------------------------------------------------------%
        function [x, w, n] = makeCubature(cubature)
            
            % Dimension of cubature
            dim    = 2;
            G      = cubature.G;
            % Basis functions used in moment-fitting
            basis  = dgBasis(dim, cubature.prescision, 'legendre');
            psi    = basis.psi;
            nDof   = basis.nDof;
            
            if cubature.prescision <= 1
                [x, ~, n] = getSquareCubaturePointsAndWeights(cubature.prescision);
                % Precision is 1, use midpoint rule
                if G.griddim == 2
                    w = G.cells.volumes;
                    x = repmat(x, G.cells.num, 1);
                    n = repmat(n, G.cells.num, 1);
                    type = 'volume';
                else
                    w = G.faces.areas;
                    type = 'face';
                    x = repmat(x, G.faces.num, 1);
                    n = repmat(n, G.faces.num, 1);
                end
            else
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % The starting point is a quadrature for a reference square
                % with more than nDof points
                if 1
                    G1 = computeGeometry(cartGrid([1,1], [2,2]));
                    G1.nodes.coords = G1.nodes.coords - 1;
                    G1 = computeVEMGeometry(G1);
                    G1 = computeCellDimensions2(G1);
                    cubTri = TriangleCubature(G1, cubature.prescision, cubature.internalConn);
                    [~, x, ~, cellNo, ~] = cubTri.getCubature(1, 'volume');
                    x = cubTri.transformCoords(x, cellNo);
                    x = unique(x, 'rows');
                else
                    k = 0;
                    nS = 0;
                    while nS < nDof
                        [x, ~, nS] = getSquareCubaturePointsAndWeights(cubature.prescision + k);
                        k = k+1;
                    end
                end
                n = size(x,1);
                % Cubature is either for faces or cells
                if G.griddim > cubature.dim
                    type = 'face';
                    elements = 1:G.faces.num;
%                     vol = G.faces.areas;
                else
                    type = 'volume';
                    elements = 1:G.cells.num;
%                     vol = G.cells.volumes;
                end

                % We use known cubature to calculate the moments
                if 0
                    if isfield(G, 'parent')
                        knownCub = CoarseGrid2DCubature(G, cubature.prescision, cubature.internalConn);
                    else
                        knownCub = TriangleCubature(G, cubature.prescision, cubature.internalConn);
                    end
                else
                    knownCub = TriangleCubature(G, cubature.prescision, cubature.internalConn);
                end
                [~, xq, wq, cellNo, faceNo] = knownCub.getCubature(elements, type);
                % Map cubature points to reference coordinates
                if G.griddim == 3
                    % Map to face reference coordinates
                    vec1  = G.faces.coordSys{1}(faceNo,:);
                    vec2  = G.faces.coordSys{2}(faceNo,:);
                    xq    = xq - G.faces.centroids(faceNo,:);
                    xq    = [sum(xq.*vec1,2), sum(xq.*vec2, 2)];
                    xq    = xq./(G.faces.dx(faceNo,:)/2);
                    count = faceNo;
                    num   = G.faces.num;
                else
                    % Map to cell reference coordiantes
                    xq    = cubature.transformCoords(xq, cellNo);
                    count = cellNo;
                    num   = G.cells.num;
                end
                % Moments
                M = cellfun(@(p) accumarray(count, wq.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, num);
                tol = eps(mean(G.cells.volumes));
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                end
                [x,w,n] = fitMoments2(x,basis, M, 'equal', G.equal);
                if numel(w) == 1
                    w = repmat(w{:}, num, 1);
                    x = repmat(x{:}, num, 1);
                    n = n*ones(num,1);
                else
                    w = vertcat(w{:});
                    x = vertcat(x{:});
                end
%                 moments = rhs(:);
%                 moments = moments./rldecode(vol, nDof, 1);
%                 [x,w,n] = fitMoments(x, basis, moments, num);
%                 w = w.*rldecode(vol, n, 1);
            end
            
            % Map from reference to physical coordinates
            if strcmp(type, 'face')
                % Face coordinates
%                 faceNo = reshape(repmat((1:G.faces.num), n, 1), [], 1);
                faceNo = rldecode((1:G.faces.num)', n, 1);
                vec1   = G.faces.coordSys{1}(faceNo,:);
                vec2   = G.faces.coordSys{2}(faceNo,:);
%                 x = x.*(G.faces.dx(faceNo,:)/2);
%                 x = repmat(x, G.faces.num, 1).*(G.faces.dx(faceNo,:)/2);
                x = x(:,1).*vec1 + x(:,2).*vec2;
                [xMin, xMax] = getMinMax(G.nodes.coords(G.faces.nodes,:), diff(G.faces.nodePos));
                dx = xMax - xMin;
                x = x.*dx(faceNo,:)/2;
                x = x + G.faces.centroids(faceNo,:);
            else
                % Cell coordinates
%                 cellNo = reshape(repmat((1:G.cells.num), n, 1), [], 1);
                cellNo = rldecode((1:G.cells.num)', n, 1);
%                 x = repmat(x, G.cells.num, 1);
                x = cubature.transformCoords(x, cellNo, true);
            end
            
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