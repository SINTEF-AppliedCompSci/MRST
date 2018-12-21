classdef CoarseGrid3DCubature < Cubature
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = CoarseGrid3DCubature(G, prescision, internalConn)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature     = cubature@Cubature(G, prescision, internalConn);
            cubature.dim = 3;
            % Make cubature points and weights
            [x, w, pos] = cubature.makeCubature();
            % Assing properties
            cubature.points  = x;
            cubature.weights = w;
            cubature.pos     = pos;
            
        end
        
        function [x, w, pos] = makeCubature(cubature)
           
            % Get cubature for parent grid
            parentCub = TetrahedronCubature(cubature.G.parent    , ...
                                            cubature.prescision  , ...
                                            cubature.internalConn);
            
            % Sort cubature points and weights according to partition
            [p, order] = sort(cubature.G.partition);
            [~, x, w]  = parentCub.getCubature(order, 'volume');
           
            np  = diff(parentCub.pos);            
            np  = accumarray(p, np(order));
            pos = [0;cumsum(np)] + 1;
                        
        end
        
    end
    
end