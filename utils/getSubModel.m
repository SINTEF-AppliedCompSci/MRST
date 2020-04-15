function [submodel, mappings] = getSubModel(model, cells, varargin)
    % Make a simulation-ready submodel defined by a subset of the cells in
    % the full model
    opt = struct('overlap', 0    , ...
                 'verbose', false);
    opt = merge_options(opt, varargin{:});
    % Get cell mappings
    cellMap = getSubCells(model, cells, opt);
    % Get face mappings
    faceMap = getSubFaces(model, cellMap);
    % Get submodel operators
    subop = getSubOperators(model, cellMap, faceMap);
    % Get subgrid
    subG = getSubGrid(model, cellMap, faceMap);
    % Make submodel
    submodel = makeSubModel(model, subG, subop, cellMap, opt);
    % Trim mappings and make mapping struct
%     faceMap  = rmfield(faceMap, 'activeConn');
    faceMap  = rmfield(faceMap, 'neighbors');
    mappings = struct('cells', cellMap, ...
                      'faces', faceMap);
end

%-------------------------------------------------------------------------%
function mappings = getSubCells(model, cells, opt)
    % Add overlap and identify external cells
    [internal, overlap] = deal(cells);
    M = getConnectivityMatrix(model.operators.N);
    for i = 1:opt.overlap
        overlap = overlap | M*overlap;
    end
    external = overlap | M*overlap;
    external = external & ~overlap;
    overlap  = overlap & ~internal;
    keep     = internal | overlap | external;
    % Construc global-to-local mapping
    renum       = zeros(model.G.cells.num,1);
    renum(keep) = 1:nnz(keep);
    % Make mapping struct
    mappings = struct('internal', internal, ...
                      'overlap' , overlap , ...
                      'external', external, ...
                      'keep'    , keep    , ...
                      'renum'   , renum   );
end

%-------------------------------------------------------------------------%
function mappings = getSubFaces(model, cellMap)
    % Find faces in internal, overlap and external regions
    type = {'internal', 'overlap', 'external'};
    N = model.G.faces.neighbors + 1;
    mappings = struct();
    for i = 1:numel(type)
        map   = [false; cellMap.(type{i})];
        faces = any(map(N),2);
        for j = 1:i-1
            faces = faces & ~mappings.(type{j});
        end
        mappings.(type{i}) = faces;
    end
    mappings.keep = mappings.internal | mappings.overlap | mappings.external;
    % Find acctive internal connections
    mappings.activeConn = all(cellMap.keep(model.operators.N), 2);
    % Construc global-to-local mapping
    renum = zeros(model.G.faces.num,1);
    renum(mappings.keep) = 1:nnz(mappings.keep);
    mappings.renum = renum;
    % Make neighbor matrix
    map = [0; cellMap.renum];
    mappings.neighbors = map(model.G.faces.neighbors(mappings.keep, :)+1);
end

%-------------------------------------------------------------------------%
function subop = getSubOperators(model, cellMap, faceMap)
    [op, subop] = deal(model.operators);
    % Dimensions
    nf = nnz(faceMap.activeConn);
    nc = nnz(cellMap.keep);
    % Transmissibility
    subop.T     = op.T(faceMap.activeConn);
    subop.T_all = op.T_all(faceMap.keep);
    % Neighboring connections
    subop.N = cellMap.renum(op.N(faceMap.activeConn, :));
    % Div and grad
    subop.M      = sparse((1:nf)'*[1, 1], subop.N, ones(nf,2)*0.5    , nf, nc);
    subop.C      = sparse((1:nf)'*[1, 1], subop.N, ones(nf,1)*[1, -1], nf, nc);
    subop.Grad   = @(x) -subop.C*x;
    subop.Div    = @(x) subop.C'*x;
    subop.AccDiv = @(acc, flux) acc + subop.C'*flux;
    if nf == 0
        % We have zero faces. Account for Matlab's preference for
        % reducing expressions of type a + [] to [].
        subop.AccDiv = @(acc, flux) acc;
        subop.Div    = @(x) zeros(nc, 1); 
    end
    if isfield(subop, 'gdz')
        subop = rmfield(subop, 'gdz');
    end
    % Face reconstructions
    subop.faceAvg   = @(x) subop.M*x;
    subop.faceUpstr = @(flag, x) faceUpstr(flag, x, subop.N, [nf, nc]);
    subop.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);
    % Pore volume
    subop.pv = op.pv(cellMap.keep);
    subop.internalConn = all(faceMap.neighbors ~= 0, 2);
    % Remove field
    if isfield(subop, 'diag_updated')
        subop = rmfield(subop, 'diag_updated');
    end
end

%-------------------------------------------------------------------------%
function subG = getSubGrid(model, cellMap, faceMap)
    % Grid
    subG        = model.G;
    subG.subset = cellMap.keep;
    % Grid cells
    subG.cells.num       = nnz(cellMap.keep);
    subG.cells.centroids = subG.cells.centroids(cellMap.keep, :);
    cells = find(cellMap.keep);
    subG.cells.faces = faceMap.renum(                             ...
               subG.cells.faces(mcolon(subG.cells.facePos(cells), ...
                                       subG.cells.facePos(cells+1)-1, 1)));
    ncf = diff(subG.cells.facePos);
    subG.cells.facePos = cumsum([0; ncf(cellMap.keep)]) + 1;
    % Grid faces
    subG.faces.num       = nnz(faceMap.keep);
    subG.faces.centroids = subG.faces.centroids(faceMap.keep, :);
    subG.faces.normals   = subG.faces.normals(faceMap.keep, :);
    subG.faces.neighbors = faceMap.neighbors;
end

%-------------------------------------------------------------------------%
function submodel = makeSubModel(model, subG, subop, cellMap, opt)
    % Make submodel
    submodel = model;
    submodel.verbose = opt.verbose;
    % Assign grid and operators
    submodel.G = subG;
    submodel.operators = subop;
    % Update operators in diagonal autodiff backend
    if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
        submodel = submodel.AutoDiffBackend.updateDiscreteOperators(submodel);
    end
    % Restrict state function groupings
    useDefault = isempty(submodel.FlowPropertyFunctions);
%     warning('Assuming default state function groupings.');
    submodel = submodel.setupStateFunctionGroupings();
    if ~useDefault
        submodel.FlowPropertyFunctions = submodel.FlowPropertyFunctions.subset(cellMap.keep);
        submodel.PVTPropertyFunctions  = submodel.PVTPropertyFunctions.subset(cellMap.keep);
        submodel.FluxDiscretization    = submodel.FluxDiscretization.subset(cellMap.keep);
        submodel.FluxDiscretization    = replaceOperators(submodel.FluxDiscretization, submodel);
    end
end

%-------------------------------------------------------------------------%
function discretization = replaceOperators(discretization, submodel)
    % Replace face upstream operator
    faceUpstr = submodel.operators.faceUpstr;
    fun       = @(sf) sf.UpwindDiscretization.setFunctionHandle(faceUpstr);
    discretization = replaceOperator(discretization, 'UpwindDiscretization', fun);
    % Replace kgrad operator (NB: Assuming simple two-point operator here!)
    tpfa    = TwoPointFluxApproximation(submodel);
    fun     = @(sf) tpfa;
    discretization = replaceOperator(discretization, 'PermeabilityGradientDiscretization', fun);
    % Replace gradient operator
    Grad    = submodel.operators.Grad;
    fun     = @(sf) Grad;
    discretization = replaceOperator(discretization, 'Grad', fun);
end

%-------------------------------------------------------------------------%
function discretization = replaceOperator(discretization, operator, fun)
    names = discretization.getNamesOfStateFunctions()';
    for name = names
        sf = discretization.getStateFunction(name{1});
        if isprop(sf, operator)
            sf.(operator) = fun(sf);
            discretization = discretization.setStateFunction(name{1}, sf);
        end
    end
end