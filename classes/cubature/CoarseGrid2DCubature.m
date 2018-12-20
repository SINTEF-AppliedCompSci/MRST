classdef CoarseGrid2DCubature < Cubature
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = CoarseGrid2DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature     = cubature@Cubature(G, prescision, internalConn);
            cubature.dim = 2;
            % Make cubature points and weights
            [x, w, pos] = cubature.makeCubature();
            % Assing properties
            cubature.points  = x;
            cubature.weights = w;
            cubature.pos     = pos;
            
        end
        
        function [x, w, pos] = makeCubature(cubature)
           
            % Get cubature for parent grid
            parentCub = TriangleCubature(cubature.G.parent    , ...
                                         cubature.prescision  , ...
                                         cubature.internalConn);
            
            % Sort cubature points and weights according to partition
            if cubature.G.griddim == 2
                part = cubature.G.partition;
                [p, order] = sort(part);
                type = 'volume';
            else
                order = cubature.G.faces.fconn;
                p     = rldecode((1:cubature.G.faces.num)', diff(cubature.G.faces.connPos), 1);
                type  = 'face';
            end
            [~, x, w]  = parentCub.getCubature(order, type);

            np  = diff(parentCub.pos);            
            np  = accumarray(p, np(order));
            pos = [0;cumsum(np)] + 1;
                        
        end
        
    end
    
end