classdef Cubature
    % Cubature class for cubatures on MRST grids
    
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
            % Setup up cubature properties. Acutual construction of
            % cubature is handled by children class.
            
            cubature.G            = G;
            cubature.prescision   = prescision;
            cubature.points       = [];
            cubature.weights      = [];
            cubature.dim          = G.griddim;
            cubature.internalConn = internalConn;
            
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature(cubature, elements, type, varargin)
            % Get cubature for given elements (either a set of cells or a
            % set of faces) of the grid.
            %
            % PARAMETERS:
            %   elements - Vector of cells or faces we want cubature for
            %   type     - String equal to either
            %              * 'volume' : Volume integral, elements
            %                           interpreted as cells
            %              * 'face'   : Surface integral, elements
            %                           interpreted as faces
            %              * 'surface': Surface integral over all faces of
            %                           elements, which are interpreted as
            %                           cells
            %
            % OPTIONAL PARAMETERS:
            %   excludeBoundary - Exclude boundary faces in cell surface
            %                     integrals
            %   internalConn    - Boolean indicating what faces are
            %                     internal connections, and therefore
            %                     included in surface integrals. If empty,
            %                     calculated from G.faces.neighbors
            %   outwardNormal   - Surface integrals over cells are often on
            %                     the form \int(v\dot n). outwardNormal =
            %                     true changes signs based on face
            %                     orientation so that normal points
            %                     outwards.
            %
            % RETURNS:
            %   W - Integration matrix - integral = W*integrand(x)
            %   x - Integration points
            %   cellNo, faceNo - cell/face number telling what cubature
            %       point goes where
            
            opt = struct('excludeBoundary', false, ...
                         'internalConn'   , []   , ...
                         'outwardNormal'  , false);         
            opt = merge_options(opt, varargin{:});
            
            if opt.excludeBoundary && isempty(opt.internalConn)
                opt.internalConn = all(cubature.G.faces.neighbors ~= 0, 2);
            end
            
            
            % Transpose if elements is row vector
            if size(elements,1) == 1, elements = elements'; end
            
            sgn = 1;
            switch cubature.G.griddim - cubature.dim
                case 0
                    % If cubature and grid dims are equal, we are dealing
                    % with a volume cubature
                    switch type
                        case 'volume'
                            % Index of cubature points and weights
                            ix = mcolon(cubature.parentPos(elements), ...
                                         cubature.parentPos(elements+1)-1);
                            % Number of cubature points for each cell
                            nq = diff(cubature.parentPos);
                            nq = nq(elements);
                            % Decode cell numbers so we know what cubature
                            % point goes where
                            cellNo = rldecode(elements, nq, 1);
                            faceNo = [];
                            
                        otherwise
                            error(['Cubature dimension tells me this', ...
                                         ' should be a volume cubature!']);
                    end
                case 1
                    % This must be a cubature for the grid faces
                    nqf = diff(cubature.parentPos); % Cub points per face
                    switch type
                        case 'face'
                            % Index of cubature points and weights
                            ix = mcolon(cubature.parentPos(elements), ...
                                         cubature.parentPos(elements+1)-1);
                            % Number of points per face
                            nq = nqf(elements);
                            % Decode face numbers so we know what cubature
                            % point goes where
                            faceNo = rldecode(elements, nq, 1);
                            cellNo = [];

                        case 'surface'
                            % Input elements are cells, we are requesting
                            % cubature for cell surfaces
                            
                            % Get all faces of cells
                            fPos = mcolon(cubature.G.cells.facePos(elements), ...
                                         cubature.G.cells.facePos(elements+1)-1);
                            faces = cubature.G.cells.faces(fPos);
                            % Number of cell faces
                            ncf = diff(cubature.G.cells.facePos);
                            ncf = ncf(elements);
                            if opt.excludeBoundary
                                % Ensure boundary faces are not counted in
                                % number of cell faces
                                cNo = rldecode((1:numel(elements))', ncf, 1);
                                ncf = ncf - accumarray(cNo, ~opt.internalConn(faces));
                                % Exclude boundary faces from faces
                                faces = faces(opt.internalConn(faces));
                            end
                            if size(faces,1) == 1, faces = faces'; end
                            % Index of cubature points and weights
                            ix = mcolon(cubature.parentPos(faces), ...
                                        cubature.parentPos(faces+1)-1);
                            % Get number of cubature points for each cell
                            nqf = diff(cubature.parentPos);
                            nq = accumarray(rldecode((1:numel(elements))', ncf, 1), nqf(faces));
                            % Avoid returning empty nq
                            if isempty(nq), nq = 0; end
                            % Decode cell and face numbers so we know what
                            % cubature point goes where
                            cellNo = rldecode(elements, nq);
                            faceNo = rldecode(faces, nqf(faces), 1);
                            if opt.outwardNormal
                                % We are computing integrals on the form
                                % \int(v\dot n), ensure correct orientation
                                sgn = 1 - 2*(cubature.G.faces.neighbors(faceNo,1) ~= cellNo);
                            end
                    end     
            end
            % Get cubature points and weights
            x = cubature.points(ix, :);
            w = cubature.weights(ix).*sgn;
            % Build integration matrix
            [ii, jj] = blockDiagIndex(ones(numel(elements),1), nq);
            W = sparse(ii, jj, w);
            
        end
            
    end
    
end