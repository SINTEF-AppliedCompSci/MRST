classdef DGDiscretization < SpatialDiscretization
    
    properties
        degree = 1          % Degree of discretization, dG(degree)
        basis  = 'legendre' % Type of basis functions. Standard is tensor 
                            % products of Legendre polynomials.
                            
        useMomentFitting    % Bool to tell the method to use experimental 
                            % unstructured cubature class
                            
        cellCubature        % Cubature for cell integrals
        faceCubature        % Cubature for face integrals
        velocityInterp      % Function for mapping face fluxes to cell
                            % velocity/ies
                            
        upwindType          % Type of upwind calculation
        internalConnParent  % If we only solve on subset of full grid, we
                            % must keep tract of internal connections in
                            % the full grid
                            
        nDof                % Number of dofs per cell
        dofPos              % Map into state.(name)dof
        sample              % AD sample
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function disc = DGDiscretization(G, varargin)
            % Constructor
            disc = disc@SpatialDiscretization(G, varargin{:});
            % Standard dG properties
            disc.degree = 1;
            disc.basis  = 'legendre';
            % Cubature
            disc.useMomentFitting = false;
            % Specifics for reordering
            disc.internalConnParent  = disc.internalConn;
            disc.G.parent            = G;
            % Merge options
            [disc, basisArgs] = merge_options(disc, varargin{:});
            % Replace basis string by dG basis object
            if ~isfield(disc.basis, 'psi')
                disc.basis = dgBasis(disc.dim, disc.degree, disc.basis, basisArgs{:});
            end
            disc.degree = disc.basis.degree;
            % Set up velocity interpolation
            disc.velocityInterp = velocityInterpolation(G, 'mimetic');
            disc.upwindType     = 'potential';
            % Create cubatures
            prescision = 2*max(disc.degree);
            isCoarse   = isfield(G, 'parent');
            if G.griddim == 2
                if isCoarse
                    volCub = CoarseGrid2DCubature(G, prescision);
                else
                    if all(disc.degree == 0) || disc.useMomentFitting
                        volCub = MomentFitting2DCubature(G, prescision);
                    else
                        volCub = TriangleCubature(G, prescision);
                    end
                end
                surfCub = LineCubature(G, prescision);
            else
                if isCoarse 
                    if ~disc.useMomentFitting
                        volCub  = CoarseGrid3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    else
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    end
                else
                    if all(disc.degree == 0)
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = MomentFitting2DCubature(G, prescision);
                    elseif disc.useMomentFitting
                        volCub  = MomentFitting3DCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    else
                        volCub  = TetrahedronCubature(G, prescision);
                        surfCub = TriangleCubature(G, prescision);
                    end
                end
            end
            disc.cellCubature  = volCub;
            disc.faceCubature = surfCub;
        end
        
        %-----------------------------------------------------------------%
        function [state, disc] = updateDofPos(disc, state)
            % Update dosfPos (position of dofs in state.sdof) based on
            % changes in state.degree, or create dofPos vector if it does
            % not exist. Dofs for cell i are found in
            %
            %   state.sdof(dofPos(:,i),:),
            %
            % Zeros are included to easily map dofs from one timestep to
            % the next.
            dp = reshape((1:disc.G.cells.num*disc.basis.nDof)', disc.basis.nDof, []);
            if isfield(state, 'nDof')
                nd = state.nDof;
            else
                nd = disc.getnDof(state);
            end
            subt = cumsum([0; disc.basis.nDof - nd(1:end-1)]);
            [ii, jj, v] = find(dp);
            % Fix dimensions
            if size(ii,1) == 1, ii = ii'; end
            if size(jj,1) == 1, jj = jj'; end
            if size(v ,1) == 1, v  =  v'; end
            % Compute position
            cnDof = cumsum(nd);
            v = v - subt(jj);
            v(v > cnDof(jj)) = 0;
            dp = full(sparse(ii, jj, v));
            % Assign
            state.nDof   = nd;
            state.dofPos = dp;
            disc.nDof    = nd;
            disc.dofPos  = dp;
        end
        
        %-----------------------------------------------------------------%
        function ix = getDofIx(disc, state, dofNo, cells, includezero)
            % Get position of dofs in state.sdof for a given cell
            %
            % PARAMETERS:
            %   state       - State with field sdof
            %   dofNo       - Dof number we want the position of. Empty
            %                 dofNo returns position of all dofs for cells
            %   cells       - Cells we want the dof position for. If empty,
            %                 positions of dofNo are returned for all cells.
            %   includeZero - Boolean indicating of we should include zeros
            %                 or not (see updateDofPos)
            %
            % RETURNS:
            %   ix - Indices into states.sdof. Dof number dofNo for cells
            %        are found in state.sdof(ix,:);
            G = disc.G;
            if nargin < 3 || (numel(dofNo) == 1 && dofNo == Inf)
                % dofNo not given, return ix for all dofs
                dofNo = 1:disc.basis.nDof;
                cells = 1:G.cells.num;
            elseif nargin < 4 || (numel(cells) == 1 && cells == Inf)
                % Cells not given, return ix for all cells
                cells = 1:G.cells.num;
            end
            ix = state.dofPos(dofNo, cells);
            ix = ix(:);
            if nargin < 5
                includezero = false;
            end
            if ~includezero
                ix(ix == 0) = [];
            end
        end
        
        %-----------------------------------------------------------------%
        function [xhat, translation, scaling] = transformCoords(disc, x, cells, inverse)
            % Transfor coordinates from physical to reference coordinates
            if nargin < 4, inverse = false; end
            [xhat, translation, scaling] = disc.cellCubature.transformCoords(x, cells, inverse);
        end
        
        %-----------------------------------------------------------------%
        function nDof = getnDof(disc, state)
            % Get number of dofs for each cell from degree.
            if nargin == 1
                deg = repmat(disc.degree, disc.G.cells.num, 1);
            else
                deg = state.degree;
            end
            if disc.degree < 0
                nDof = 0;
            else
                nDof = factorial(deg + disc.dim)...
                       ./(factorial(disc.dim).*factorial(deg));
            end
        end
        
        %-----------------------------------------------------------------%
        function state = mapDofs(disc, state, state0, name)
            % Map dofs from state0 to state, typically from one timestep to
            % the next, when we start with maximum number of dofs in all
            % cells.
            state = disc.updateDofPos(state);
            if nargin == 3
                name = 'sdof';
            else
                name = [name, 'dof'];
            end
            if all(state.nDof == state0.nDof)
                state.(name) = state0.(name);
            else
                dof = zeros(sum(state.nDof), size(state0.(name),2));
                for dofNo = 1:disc.basis.nDof
                    % We may be solving only a subset of the gridcells, so
                    % we include zeros in dofIx to keep track of where old
                    % dofs maps to new ones.
                    ix  = disc.getDofIx(state , dofNo, (1:disc.G.cells.num)', true);
                    ix0 = disc.getDofIx(state0, dofNo, (1:disc.G.cells.num)', true);
                    dof(ix(ix0 > 0 & ix > 0),:) = state0.(name)(ix0(ix0 > 0 & ix > 0),:);
                end
                state.(name) = dof;
            end
        end
        
        %-----------------------------------------------------------------%
        function p = evaluateDGVariable(disc, x, cells, state, pdof, psi)
            % Evaluate dG variable at cubature points 'x' in 'cells'
            if isempty(pdof)
                return
            else
                if nargin < 6 || isempty(psi)
                    x    = disc.transformCoords(x, cells);
                    psi  = disc.basis.psi;
                    getx = @(keep) x(keep,:);
                else
                    getx = @(keep) keep;
                end
                p = pdof(cells,:)*0;
                for dofNo = 1:disc.basis.nDof
                    keep = state.nDof(cells) >= dofNo;
                    ix = disc.getDofIx(state, dofNo, cells(keep));
                    if all(keep)
                        p = p + pdof(ix,:).*psi{dofNo}(getx(keep));
                    else
                        p(keep, :) = p(keep, :) + pdof(ix,:).*psi{dofNo}(getx(keep));
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function p = evaluateProp(disc, state, dof, type, elements)
            % Convenience function for evaluation of model properties.
            % Called by e.g., TransportModelDG.
            psi = [];
            if nargin < 5
                elements = Inf;
            end
            switch type
                case 'cell'
                    if isfield(state, 'psi_c')
                        psi = state.psi_c;
                        cNo = state.cells;
                        x   = [];
                    else
                        [~ , x, cNo] = disc.getCubature(elements, 'cell');
                    end
                case 'face'
                    if isfield(state, 'psi_f')
                        psi = state.psi_f;
                        cNo = state.fcells;
                        x   = [];
                    else
                        [~ , x, ~, fNo] = disc.getCubature(elements, 'face');
                        x   = repmat(x, 2, 1);
                        N   = disc.G.faces.neighbors;
                        cNo = [N(fNo,1); N(fNo,2)];
                    end
            end
            if ~isempty(psi) && ~iscell(psi)
                psi = mat2cell(psi, size(psi,1), ones(1,disc.basis.nDof));
            end
            if iscell(dof)
                p = cell(1,numel(dof));
                for i = 1:numel(dof)
                    p{i} = disc.evaluateDGVariable(x, cNo, state, dof{i}, psi);
                end
            else
                p = disc.evaluateDGVariable(x, cNo, state, dof, psi);
            end
        end
        
        %-----------------------------------------------------------------%
        function state = evaluateBasisFunctions(disc, state, cells, faces)
            % Evaluate basis functions in cubature points
            [~ , xc, cNo     ] = disc.getCubature(cells, 'cell');
            [~ , xf, ~  , fNo] = disc.getCubature(faces, 'face');
            N    = disc.G.faces.neighbors;
            if all(any(N(fNo,:) == 0,2))
                fcNo = sum(N(fNo,:),2);
            else
                xf   = repmat(xf, 2, 1);
                fcNo = [N(fNo,1); N(fNo,2)];
            end
            xc = disc.transformCoords(xc, cNo );
            xf = disc.transformCoords(xf, fcNo);
            [psi_c, psi_f] = deal(disc.basis.psi');
            for dofNo = 1:disc.basis.nDof
                psi_c{dofNo} = psi_c{dofNo}(xc);
                psi_f{dofNo} = psi_f{dofNo}(xf);
            end
            state.psi_c  = psi_c;
            state.psi_f  = psi_f;
            state.cells  = cNo;
            state.fcells = fcNo;
            state.faces  = fNo;
            state.mcells = cells;
            state.mfaces = faces;
        end
        
        %-----------------------------------------------------------------%
        function fill = getFillSat(disc, state)
            % Get total saturation dofs of all cells
            fill = zeros(sum(state.nDof),1);
            ix   = disc.getDofIx(state, 1);
            fill(ix) = 1;
        end
        
        %-----------------------------------------------------------------%
        function varargout = getCellMean(disc, state, cells, varargin)
            % Get average cell value from dofs
            W = disc.getCubature(cells, 'cell');
            val = cell(numel(varargin),1);
            for i = 1:nargin-3
                dof = varargin{i};
                v = disc.evaluateProp(state, dof, 'cell', cells);
                if iscell(dof)
                    for j = 1:numel(dof)
                        v{j} = W*v{j};
                    end
                else
                    v = W*v;
                end
                val{i} = v;
            end
            varargout = val;
        end
        
        %-----------------------------------------------------------------%
        function ip = inner(disc, u, v, differential, elements)
            % Compute dG inner products
            %
            % PARAMETERS:
            %   u            - Either scalar vector or SpatialVector
            %   v            - Basis function cell array
            %   differential - 'dV': cell integral, 'dS' face integral
            %   elements     - Elements we compute innter product over
            %
            % RETURNS:
            %   ip - Inner product (u,v) over elements
            if nargin < 5
                elements = Inf;
            end 
            switch differential
                case 'dV'
                    ip = disc.cellInt(u, v, elements);
                case 'dS'
                    ip = disc.faceInt(u, v, elements);
            end
        end
        
        %-----------------------------------------------------------------%
        function I = cellInt(disc, u, v, cells)
            % Integrate integrand over cells
            nDofMax  = numel(v);     % Maximum number of dofs
            % Empty cells means all cells in grid
            if isempty(cells)
                cells = (1:disc.G.cells.num)';
            end
            % Get cubature for all cells, transform coordinates to ref space
            [W, x, cellNo, ~] = disc.getCubature(cells, 'cell');
            if isinf(cells)
                cells = (1:disc.G.cells.num)';
            end
            [x, ~, scaling]   = disc.transformCoords(x, cellNo);
            % Evaluate integrals
            I = disc.sample*0;
            for dofNo = 1:nDofMax
                keepCells = disc.nDof(cells) >= dofNo;
                if any(keepCells)
                    ix = disc.getDofIx(disc, dofNo, cells(keepCells));
                    if isa(u, 'SpatialVector')
                        i = W*dot(u,v{dofNo}(x).*scaling);
                    else
                        i = W*(u.*v{dofNo}(x));
                    end
                    I(ix) = i(keepCells);
                elseif numel(cells) == disc.G.cells.num
                    warning('No cells with %d dofs', dofNo);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function I = faceInt(disc, u, v, faces)
            % Integrate integrand over faces
            nDofMax = disc.basis.nDof; % maximum number of dofs
            % Get cubature
            [W, x, ~, faceNo] = disc.getCubature(faces, 'face');
            if isinf(faces)
                % faces = inf means all internal faces in the grid
                faces = find(disc.internalConn)';
            end
            N = disc.G.faces.neighbors;
            if isempty(faceNo)
                I = 0;
                return
            end
            I = disc.sample*0;
            facesSide = faces;
            xSide     = x;
            uSide     = u;
            WSide     = W;
            for side = 1:2
                cells  = N(faces , side);
                cellNo = N(faceNo, side);
                ix_c   = cells > 0;
                if ~all(ix_c)
                    cells     = cells(ix_c);
                    facesSide = faces(ix_c);
                    ix_cNo    = cellNo > 0;
                    cellNo    = cellNo(ix_cNo);
                    xSide     = x(ix_cNo,:);
                    uSide     = u(ix_cNo);
                    WSide     = W(ix_c,ix_cNo);
                end
                f2c = sparse(cells, (1:numel(facesSide))', 1, disc.G.cells.num, numel(facesSide));
                [xf, ~, ~] = disc.transformCoords(xSide, cellNo);
                % Evaluate integrals
                sgn = (-1).^(side-1);
                for dofNo = 1:nDofMax
                    keepCells = disc.nDof(cells) >= dofNo;
                    kc = f2c*keepCells > 0;
                    if any(keepCells)
                        ix = disc.getDofIx(disc, dofNo, kc');
                        i  = f2c*(sgn.*WSide*(uSide.*v{dofNo}(xf)));
                        I(ix) = I(ix) + i(kc);
                    elseif numel(facesSide) == nnz(disc.internalConn)
                        warning('No cells with %d dofs', dofNo);
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function [W, x, cellNo, faceNo] = getCubature(disc, elements, type, varargin)
            % Get cubature for elements. Wrapper for cubature class
            % function getCubature, with mapping of elements before and
            % after in case we are solving on a subgrid
            useMap = isfield(disc.G, 'mappings');
            if useMap
                % Map elements to old numbering
                maps = disc.G.mappings; 
                switch type 
                    case {'cell', 'surface'}
                        elements = maps.cellMap.new2old(elements);
                    case 'face'
                        elements = maps.faceMap.new2old(elements);
                end
            end
            % Get correct cubature type
            switch type 
                case 'cell'
                    cubature = disc.cellCubature; 
                case {'surface', 'face'}
                    cubature = disc.faceCubature;
            end
            % Make options struct
            opt = struct('excludeBoundary', true                   , ...
                         'internalConn'   , disc.internalConnParent, ...
                         'outwardNormal'  , true                   );
            opt = merge_options(opt, varargin{:});
            % Get cubature from cubature class
            [W, x, ~, cellNo, faceNo] = cubature.getCubature(elements, type, ...
                                 'excludeBoundary', opt.excludeBoundary    , ...
                                 'internalConn'   , opt.internalConn       , ...
                                 'outwardNormal'  , opt.outwardNormal      );
            if useMap
                % Map elements back to new numbering
                cellNo = maps.cellMap.old2new(cellNo);
                faceNo = maps.faceMap.old2new(faceNo);
            end
        end
        
        %-----------------------------------------------------------------%
        function [vMin, vMax] = getMinMax(disc, state, dof)
            % Get minimum and maximum value of dG variable for each cell.
            % Evaluated in cubature points
            G = disc.G;
            % Get all quadrature points for all cells
            [~, xF, ~, fF] = disc.getCubature((1:G.faces.num), 'face');
            xF = repmat(xF, 2, 1);
            cF = [G.faces.neighbors(fF,1); G.faces.neighbors(fF,2)];
            ix = cF == 0;
            cF(ix) = [];
            xF(ix,:) = [];
            [cF, ix] = sort(cF);
            xF = xF(ix,:);
            [~, xC, cC, ~] = disc.getCubature(Inf, 'cell');
            % Evaluate saturation at faces
            vF = disc.evaluateDGVariable(xF, cF, state, dof);
            [~, nF] = rlencode(cF);
            [vMinF, vMaxF] = getMinMax(vF, nF);
            % Evaluate saturation at cells
            vC = disc.evaluateDGVariable(xC, cC, state, dof);
            [~, nC] = rlencode(cC);
            % Find min/max saturation
            [vMinC, vMaxC] = getMinMax(vC, nC);
            vMin = min(vMinF, vMinC);
            vMax = max(vMaxF, vMaxC);
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
