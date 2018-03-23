classdef Cubature
    
    properties
        
        G
        prescision
        points
        weights
        numPoints
        dim
        internalConn
        
    end
    
    methods
        
        function cub = Cubature(G, prescision, internalConn)
            
            cub.G = G;
            cub.prescision = prescision;
            cub.points  = [];
            cub.weights = [];
            cub.dim = G.griddim;
            cub.internalConn = internalConn;
            
        end

        function [x, w, ii, jj, cellNo, faceNo] = getCubature(cub, elements, type)
            
            if size(elements,1) == 1, elements = elements'; end
            g = cub.G;
            
            switch g.griddim - cub.dim
                case 0
                    
                    ix = mcolon(cub.parentPos(elements), cub.parentPos(elements+1)-1);
                    
                    nq = diff(cub.parentPos);
                    nq = nq(elements);
                    
                    cellNo = rldecode(elements, nq, 1);
                    faceNo = [];
                    sign   = 1;
                    
                case 1
                    
                    if strcmp(type, 'face')
                        faces = (1:g.faces.num)';
                        elements = faces;
                    else
                        faces = g.cells.faces(mcolon(g.cells.facePos(elements), g.cells.facePos(elements+1)-1));
                    end
                    
                    if size(faces,1) == 1, faces = faces'; end
                    ix = mcolon(cub.parentPos(faces), cub.parentPos(faces+1)-1);
                    
                    nqf = diff(cub.parentPos);
                    faceNo = rldecode(faces, nqf(faces), 1);
                    if strcmp(type, 'face')
                        nq = nqf;
                        cellNo = [];
                        sign = 1;
                    else
                        faces = faces(cub.internalConn(faces));
                        ncf = accumarray(rldecode((1:g.cells.num)', diff(g.cells.facePos), 1), cub.internalConn(faces));
                        nq = accumarray(rldecode(elements, ncf(elements), 1), nqf(faces));
                        cellNo = rldecode(elements, nq);
                        sign   = 1 - 2*(g.faces.neighbors(faceNo,1) ~= cellNo);
                    end
                    
            end
            
            x = cub.points(ix, :);
            w = cub.weights(ix).*sign;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
                        
        end
            
    end
    
end