classdef Cubature
    
    properties
        
        G
        prescision
        points
        weights
        numPoints
        dim
        internalConn
        W
        
    end
    
    methods
        
        function cub = Cubature(G, prescision, internalConn)
            
            cub.G = G;
            cub.prescision = prescision;
            cub.points  = [];
            cub.weights = [];
            cub.dim = G.griddim;
            cub.internalConn = internalConn;
            W = [];
            
        end

        function [W, x, w, ii, jj, cellNo, faceNo] = getCubature(cub, elements, type)
            
            if size(elements,1) == 1, elements = elements'; end
            g = cub.G;
            
            switch g.griddim - cub.dim
                case 0
                    
                    ix = mcolon(cub.parentPos(elements), cub.parentPos(elements+1)-1);
                    
                    nq = diff(cub.parentPos);
                    nq = nq(elements);
                    
                    cellNo = rldecode(elements, nq, 1);
                    faceNo = [];
                    sgn   = 1;
                    
                case 1
                    
                    
                    switch type
                        case 'face'
                            faces = elements;
%                             faces = (1:g.faces.num)';
%                             elements = faces;
                            ix = mcolon(cub.parentPos(faces), cub.parentPos(faces+1)-1);
                        case 'cellsurface'
                            faces = g.cells.faces(mcolon(g.cells.facePos(elements), g.cells.facePos(elements+1)-1));
                    end
                    
                    if size(faces,1) == 1, faces = faces'; end
                    
                    nqf = diff(cub.parentPos);
                    switch type
                        case 'face'
                            faceNo = rldecode(faces, nqf(faces), 1);
                            nq = nqf(faces);
                            cellNo = [];
                            sgn = 1;
                        case 'cellsurface'
                            ncf = accumarray(rldecode((1:g.cells.num)', diff(g.cells.facePos), 1), cub.internalConn(faces));
                            faces = faces(cub.internalConn(faces));
                            ix = mcolon(cub.parentPos(faces), cub.parentPos(faces+1)-1);
                            nq = accumarray(rldecode(elements, ncf(elements), 1), nqf(faces));
                            faceNo = rldecode(faces, nqf(faces), 1);
                            cellNo = rldecode(elements, nq);
                            sgn   = 1 - 2*(g.faces.neighbors(faceNo,1) ~= cellNo);
                    end
                    
            end
            
            x = cub.points(ix, :);
            w = cub.weights(ix).*sgn;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
            
        end
            
    end
    
end