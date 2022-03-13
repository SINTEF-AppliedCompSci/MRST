classdef Cubature
    % Cubature class for cubatures on MRST grids
    
    properties
        G            % Grid the cubature is defined on
        prescision   % Prescision of cubature
        points       % Cubature points
        weights      % Cubature weights
        numPoints    % Number of cubature points
        pos          % Cubature points/weights for element e can be found
                     % in ix = cubature.pos(e):cubature.pos(e+1)-1
        dim          % Dimension of cubature
        internalConn % Boolean indicating if face is internal connection
    end
    
    properties (Access = protected)
        fullCubature
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = Cubature(G, prescision, dim, varargin)
            % Setup up cubature properties. Acutual construction of
            % cubature is handled by children class.
            cubature.G          = G;
            cubature.prescision = prescision;
            cubature.dim          = dim;
            internalConn = all(G.faces.neighbors > 0, 2);
            cubature.internalConn = internalConn;
            cubature = merge_options(cubature, varargin{:});
            % Make cubature points and weights
            cubature = cubature.makeCubature();
            % Get cubature for entire domain for quick access
            switch G.griddim - cubature.dim
                case 0
                    elements = (1:G.cells.num)';
                    type     = 'cell';
                case 1
                    elements = find(cubature.internalConn);
                    type     = 'face';
            end
            [W, x, w, cellNo, faceNo] = cubature.getCubature(elements, type);
            cubature.fullCubature = struct('W', W, 'x', x, 'w', w, ...
                                       'cellNo', cellNo, 'faceNo', faceNo);
        end
        
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature) %#ok
            error('Base class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function cubature = assignCubaturePointsAndWeights(cubature, x, w, n)
            % Assign points and weights
            switch cubature.G.griddim - cubature.dim
                case 0
                    w = w./rldecode(cubature.G.cells.volumes, n, 1);
                case 1
                    w = w./rldecode(cubature.G.faces.areas, n, 1);
            end
            cubature.points    = x;
            cubature.weights   = w;
            cubature.numPoints = n;
            % Construct cubature position vector
            cubature.pos = [0; cumsum(n)] + 1;
        end
            
        %-----------------------------------------------------------------%
        function [W, x, w, cellNo, faceNo] = getCubature(cubature, elements, type, varargin)
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
            
            if isinf(elements)
                W      = cubature.fullCubature.W;
                x      = cubature.fullCubature.x;
                w      = cubature.fullCubature.w;
                cellNo = cubature.fullCubature.cellNo;
                faceNo = cubature.fullCubature.faceNo;
                return
            end
            
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
                        case 'cell'
                            % Index of cubature points and weights
                            ix = mcolon(cubature.pos(elements), ...
                                        cubature.pos(elements+1)-1);
                            % Number of cubature points for each cell
                            nq = diff(cubature.pos);
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
                    nqf = diff(cubature.pos); % Cub points per face
                    switch type
                        case 'face'
                            % Index of cubature points and weights
                            ix = mcolon(cubature.pos(elements), ...
                                        cubature.pos(elements+1)-1);
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
                            ix = mcolon(cubature.pos(faces), ...
                                        cubature.pos(faces+1)-1);
                            % Get number of cubature points for each cell
                            nqf = diff(cubature.pos);
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
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(cub, x, cells, inverse)
            % Transfor coordinates from physical to reference coordinates. 
            %
            % PARAMETERS:
            %   x         - Coordinates in physical space
            %   cells     - Cells we want reference coordinates for, cells(ix)
            %               are used when transforming x(ix,:)
            %   inverse   - Boolean indicatiing if we are mapping 
            %               to (inverse = false) or from (inverse = true)
            %               reference coordiantes. Default = false.
            %   useParent - Boolean indicating if we are working on the
            %               full grid (G.parent) or a subgrid.
            %
            % RETURNS:
            %   xhat        - Transformed coordinates
            %   translation - Translation applied to x
            %   scaling     - Scaling applied to x
            
            if nargin < 4, inverse = false; end
            
            % Coordinates are centered in cell center
            translation = -cub.G.cells.centroids(cells,:);
            if isfield(cub.G.cells, 'dx')
                % Scaling found from dimensions of minimum bounding box
                % aligned with coordinate axes that contains the cell
                scaling = 1./(cub.G.cells.dx(cells,:)/2);
            else
                % If G.cells.dx is not computed, we use approximation
                dx = cub.G.cells.volumes(cells).^(1/cub.G.griddim);
                scaling = repmat(1./(dx/2), 1, cum.dim);
            end
            
            if ~inverse
                xhat = (x + translation).*scaling;
                xhat = xhat(:, 1:cub.dim);
                scaling     = scaling(:, 1:cub.dim);
                translation = translation(:, 1:cub.dim);
%                 assert(all(all(abs(xhat)<=1)))
            else
                xhat = x./scaling - translation;
            end
               
        end
            
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
