classdef CoarseGrid1DCubature < Cubature
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = CoarseGrid1DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, internalConn);
            cubature.dim       = 1;
            % Make cubature points and weights
            [x, w, pos] = cubature.makeCubature();
            % Assing properties
            cubature.points    = x;
            cubature.weights   = w;
            cubature.pos       = pos;
            
        end
        
        function [x, w, pos] = makeCubature(cubature)
           
            % Get cubature for parent grid
            parentCub = LineCubature(cubature.G.parent    , ...
                                     cubature.prescision  , ...
                                     cubature.internalConn);
            
            % Sort cubature points and weights according to partition
            
            order = cubature.G.faces.fconn;
            p     = rldecode((1:cubature.G.faces.num)', diff(cubature.G.faces.connPos), 1);
            [~, x, w]  = parentCub.getCubature(order, 'face');

            np  = diff(parentCub.pos);            
            np  = accumarray(p, np(order));
            pos = [0;cumsum(np)] + 1;
                        
        end
        
    end
    
end