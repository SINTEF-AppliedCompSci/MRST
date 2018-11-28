classdef LineCubature < Cubature
    % 1D line cubatureature for faces in 2D MRST grids
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = LineCubature(G, prescision, internalConn)
            % Set up line cubatureature
            
            % Make sure we have a 2D grid structure
            assert(G.griddim == 2, 'LineCubature only supported for 2D grids')
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            % Make cubatureature points and weights
            [x, w, n] = cubature.makeCubature();
            % Assign properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            cubature.dim       = 1;
            % Construct cubatureature position vector
            cubature.pos       = (0:cubature.numPoints:cubature.numPoints*G.faces.num)' + 1;
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
%             R     = zeros(nPts, G.griddim*nLin);
%             for dNo = 1:G.griddim
%                 R(:,dNo:G.griddim:end) = reshape(xn(:,dNo), nPts, []);
%             end
%             % Map coordinates
%             xx = x*R;
%             x = zeros(nLin*nq, G.griddim);
%             for dNo = 1:G.griddim
%                 xtmp     = xx(:, dNo:G.griddim:end);
%                 x(:,dNo) = xtmp(:);
%             end
            
        end
        
        %-----------------------------------------------------------------%
        function [x, w, n, linNo] = makeCubature(cubature)
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
            
        end
        
    end
    
end