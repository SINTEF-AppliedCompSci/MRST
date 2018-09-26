classdef Cubature
    
    properties
        
        G            % Grid the cubature is defined on
        prescision   % Prescision of cubature
        points       % Cubature points
        weights      % Cubature weights
        numPoints    % Number of cubature points
        dim          % Dimension of cubature
        internalConn % Boolean indicating if face is internal connection
        
    end
    
    methods
        
        function cubature = Cubature(G, prescision, internalConn)
            
            cubature.G            = G;
            cubature.prescision   = prescision;
            cubature.points       = [];
            cubature.weights      = [];
            cubature.dim          = G.griddim;
            cubature.internalConn = internalConn;
            
        end

        function [W, x, w, ii, jj, cellNo, faceNo] = getCubature(cubature, elements, type)
            
            if size(elements,1) == 1, elements = elements'; end
            g = cubature.G;
            
            switch g.griddim - cubature.dim
                case 0
                    
                    ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
                    
                    nq = diff(cubature.parentPos);
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
                            ix = mcolon(cubature.parentPos(faces), cubature.parentPos(faces+1)-1);
                        case 'cellsurface'
                            faces = g.cells.faces(mcolon(g.cells.facePos(elements), g.cells.facePos(elements+1)-1));
                    end
                    
                    if size(faces,1) == 1, faces = faces'; end
                    
                    nqf = diff(cubature.parentPos);
                    switch type
                        case 'face'
                            faceNo = rldecode(faces, nqf(faces), 1);
                            nq = nqf(faces);
                            cellNo = [];
                            sgn = 1;
                        case 'cellsurface'
                            ncf = accumarray(rldecode((1:g.cells.num)', diff(g.cells.facePos), 1), cubature.internalConn(faces));
                            faces = faces(cubature.internalConn(faces));
                            ix = mcolon(cubature.parentPos(faces), cubature.parentPos(faces+1)-1);
                            nq = accumarray(rldecode(elements, ncf(elements), 1), nqf(faces));
                            faceNo = rldecode(faces, nqf(faces), 1);
                            cellNo = rldecode(elements, nq);
                            sgn   = 1 - 2*(g.faces.neighbors(faceNo,1) ~= cellNo);
                    end
                    
            end
            
            x = cubature.points(ix, :);
            w = cubature.weights(ix).*sgn;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
            
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature2(cubature, elements, type, excludeBoundary)
            % TODO: Move mapping of elements to DGDiscretizaiton.
            if size(elements,1) == 1, elements = elements'; end

            useMap = false;
            if isfield(cubature.G, 'mappings')
                useMap = true;
                maps = cubature.G.mappings;
                cubature.G = cubature.G.parent;
                switch type 
                    case {'volume', 'surface'}
                        elements = maps.cellMap.new2old(elements);
                    case 'face'
                        elements = maps.faceMap.new2old(elements);
                end
            end
            
            switch cubature.G.griddim - cubature.dim
                case 0
                    switch type
                        case 'volume'
                            ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
                            nq = diff(cubature.parentPos);
                            nq = nq(elements);

                            cellNo = rldecode(elements, nq, 1);
                            faceNo = [];

                            sgn   = ones(numel(ix),1);
                        otherwise
                            error(['Cubature dimension is not consistent ' ...
                                   'with anything else than volume cubature']);
                    end
                case 1
                    
                    switch type
                        case 'face'
                    
                            nqf = diff(cubature.parentPos);
                            faceNo = rldecode(elements, nqf(elements), 1);
                            ix = mcolon(cubature.parentPos(elements), cubature.parentPos(elements+1)-1);
                            nq = nqf(elements);
                            cellNo = [];
                            sgn = 1;

                        case 'surface'
                    
                            faces = cubature.G.cells.faces(mcolon(cubature.G.cells.facePos(elements), cubature.G.cells.facePos(elements+1)-1));
                            ncf = accumarray(rldecode((1:cubature.G.cells.num)', diff(cubature.G.cells.facePos), 1), cubature.internalConn(cubature.G.cells.faces(:,1)));
                            if excludeBoundary
                                faces = faces(cubature.internalConn(faces));
                            end
                            if size(faces,1) == 1, faces = faces'; end
                            ix = mcolon(cubature.parentPos(faces), cubature.parentPos(faces+1)-1);

                            nqf = diff(cubature.parentPos);
                            nq = accumarray(rldecode((1:numel(elements))', ncf(elements), 1), nqf(faces));
                            if isempty(nq), nq = 0; end

                            cellNo = rldecode(elements, nq);
                            faceNo = rldecode(faces, nqf(faces), 1);
                            sgn   = 1 - 2*(cubature.G.faces.neighbors(faceNo,1) ~= cellNo);
                    end
%                     ixf = nan(numel(ix),1);
%                     ixf(disc.internalConn(faceNo)) = 1:nnz(disc.internalConn(faceNo));
                    
            end
            
            if useMap
                cellNo = maps.cellMap.old2new(cellNo);
                faceNo = maps.faceMap.old2new(faceNo);
            end
            
            x = cubature.points(ix, :);
            w = cubature.weights(ix).*sgn;
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
%             W = cubature.W(cells, ixf);
            
        end
            
    end
    
end