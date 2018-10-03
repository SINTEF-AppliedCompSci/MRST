classdef LineCubature < Cubature
    
    properties
        
        linNo
        parentPos
        
    end
    
    methods
        
        function cub = LineCubature(G, prescision, internalConn)
            
            cub = cub@Cubature(G, prescision, internalConn);
            [x, w, n, linNo] = cub.makeCubature();
            cub.points = x;
            cub.weights = w;
            cub.numPoints = n;
            cub.linNo = linNo;
            cub.dim = 1;
            
            
            cub.parentPos = (0:cub.numPoints:cub.numPoints*G.faces.num)' + 1;
            
        end
        
        function [x, w, n] = getCubaturePointsAndWeights(cub)
            
            presc = cub.prescision;
            presc = presc + (presc == 0);
            l     = legendrePolynomials(presc); 
            
            dl = cellfun(@(l) dx(l,1), l, 'unif', false);
            dl = dl{end};
            
            if presc <= 1
                x = 0;
                
            elseif presc == 2
                x = [-sqrt(1/3); sqrt(1/3)];
                
            elseif presc == 3
                x = [-sqrt(3/5); 0; sqrt(3/5)];
                
            elseif presc == 4
                x = [-sqrt(3/7 + 2/7*sqrt(6/5)); -sqrt(3/7 - 2/7*sqrt(6/5)); ...
                      sqrt(3/7 - 2/7*sqrt(6/5));  sqrt(3/7 + 2/7*sqrt(6/5))];
                
            elseif presc == 5
                a = 5; b = 2*sqrt(10/7);
                x = [-1/3*sqrt(a + b); -1/3*sqrt(a - b); 0; ...
                     +1/3*sqrt(a - b); +1/3*sqrt(a + b)];
                 
            elseif presc <= 7
                
                presc = 7;
                presc = presc + (presc == 0);
                l  = legendrePolynomials(presc);
                dl = cellfun(@(l) dx(l,1), l, 'unif', false);
                dl = dl{end};
                
                x = [-0.949107912342758524526189684048;
                     -0.741531185599394439863864773281;
                     -0.405845151377397166906606412077;
                      0;
                      0.405845151377397166906606412077;
                      0.741531185599394439863864773281;
                      0.949107912342758524526189684048];
                  
            end
            
            w = 2./((1 - x.^2).*dl(x).^2);
            
            x = (x + 1)/2;
            x = [x, 1-x];
            w = w/2;
            
            n = numel(w);
             
        end
        
        function x = mapCoords(cub, x)
            
            G = cub.G;
            npts = 2;
            nq = size(x,1);
            nodes = G.faces.nodes(:,1);
            xl = G.nodes.coords(nodes,:);
            
            nLin = G.faces.num;
            R = zeros(npts, G.griddim*nLin);
            for dNo = 1:G.griddim
                R(:,dNo:G.griddim:end) = reshape(xl(:,dNo), npts, []);
            end
            
            xx = x*R;
            x = zeros(nLin*nq, G.griddim);
            for dNo = 1:G.griddim
                xtmp = xx(:, dNo:G.griddim:end);
                x(:,dNo) = xtmp(:);
            end
            
            
        end
        
        function [x, w, n, linNo] = makeCubature(cub)
            
            G = cub.G;
            nLin = G.faces.num;
            [x, w, n] = cub.getCubaturePointsAndWeights();
            x = cub.mapCoords(x);
            
            linNo = reshape(repmat(1:G.faces.num, n, 1), [], 1);
            w = repmat(w, nLin, 1).*G.faces.areas(linNo);
            
        end
        
    end
    
end