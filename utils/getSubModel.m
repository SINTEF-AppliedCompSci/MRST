function [submodel, mappings] = getSubModel(model, cells, varargin)
    % Make a simulation-ready submodel defined by a subset of the cells in
    % the full model

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('overlap'  , 0    , ...
                 'verbose'  , false, ...
                 'plottable', true );
    opt = merge_options(opt, varargin{:});
    % Get cell mappings
    cellMap = getSubCells(model, cells, opt);
    % Get face mappings
    faceMap = getSubFaces(model, cellMap);
    % Get submodel operators
    subop = getSubOperators(model, cellMap, faceMap);
    % Get subgrid
    subG = getSubGrid(model, cellMap, faceMap, opt);
    % Get subrock
    subRock = getSubRock(model, cellMap);
    % Make submodel
    submodel = makeSubModel(model, subG, subRock, subop, cellMap, faceMap, opt);
    % Trim mappings and make mapping struct
    faceMap  = rmfield(faceMap, 'activeConn');
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
    % Transpose in case we have a single internal connection
    if size(subop.N,2) == 1, subop.N = subop.N'; end
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
function subG = getSubGrid(model, cellMap, faceMap, opt)
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
    % Grid faces if present
    if ~isfield(subG.faces, 'centroids')
        return;
    end
    subG.faces.num       = nnz(faceMap.keep);
    subG.faces.centroids = subG.faces.centroids(faceMap.keep, :);
    subG.faces.normals   = subG.faces.normals(faceMap.keep, :);
    subG.faces.areas     = subG.faces.areas(faceMap.keep);
    if isfield(subG.faces, 'tag')
        subG.faces.tag = subG.faces.tag(faceMap.keep);
    end
    subG.faces.neighbors = faceMap.neighbors;
    nfn = diff(subG.faces.nodePos);
    if opt.plottable
        faces = find(faceMap.keep);
        subG.faces.nodes = subG.faces.nodes(mcolon(subG.faces.nodePos(faces), ...
                                                   subG.faces.nodePos(faces+1)-1,1));
    end
    subG.faces.nodePos = cumsum([0; nfn(faceMap.keep)]) + 1;
end

%-------------------------------------------------------------------------%
function subRock = getSubRock(model, cellMap)
    % Extract submodel rock properties
    subRock = getSubRockFields(model.rock, cellMap);
end

%-------------------------------------------------------------------------%
function subRock = getSubRockFields(rock, cellMap)
    fnames = reshape(fieldnames(rock), 1, []);
    subRock = rock;
    for f = fnames
        if size(rock.(f{1}),1) == numel(cellMap.keep)
            subRock.(f{1}) = rock.(f{1})(cellMap.keep,:);
        elseif isstruct(rock.(f{1}))
            subRock.(f{1}) = getSubRockFields(rock.(f{1}), cellMap);
        end
    end
end

%-------------------------------------------------------------------------%
function submodel = makeSubModel(model, subG, subRock, subop, cellMap, faceMap, opt)
    % Make submodel
    submodel = model;
    submodel.verbose = opt.verbose;
    % Assign grid, rock, and operators
    submodel.G         = subG;
    submodel.rock      = subRock;
    submodel.operators = subop;
    % Update operators in diagonal autodiff backend
    if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
        submodel = submodel.AutoDiffBackend.updateDiscreteOperators(submodel);
    end
    % Restrict state function groupings if already set
    if ~isempty(submodel.getStateFunctionGroupings())
        warning('Assuming default state function groupings.');
        
%         submodel.FlowPropertyFunctions = submodel.FlowPropertyFunctions.subset(cellMap.keep);
%         submodel.FlowPropertyFunctions = submodel.FlowPropertyFunctions.subset(':');
%         if isprop(submodel.FlowPropertyFunctions.RelativePermeability, 'scalers')
%             submodel.FlowPropertyFunctions.RelativePermeability.scalers = ...
%                 getScalersSubset(submodel.FlowPropertyFunctions.RelativePermeability.scalers, cellMap.keep);
%         end
        
        % Replace flow property functions
        submodel.FlowPropertyFunctions ...
            = replaceFlowProps(submodel.FlowPropertyFunctions, cellMap);
        % Replace PVT property functions
        submodel.PVTPropertyFunctions ...
            = replacePVTProps(submodel.PVTPropertyFunctions, cellMap);
        % Replace Flow discretization
        submodel.FlowDiscretization ...
            = replaceFlowDisc(submodel.FlowDiscretization, cellMap, faceMap, submodel);
%         
%         submodel.PVTPropertyFunctions  = submodel.PVTPropertyFunctions.subset(cellMap.keep);
%         submodel.FlowDiscretization    = submodel.FlowDiscretization.subset(cellMap.keep);
%         submodel.FlowDiscretization    = replaceOperators(submodel.FlowDiscretization, submodel);
    end
end

%-------------------------------------------------------------------------%
function props = replaceFlowProps(props, map)
    props = props.subset(map.keep);
    props = props.subset(':');
    % Update relperm-specific parameters
    if isprop(props.RelativePermeability, 'scalers') ...
            && ~isempty(props.RelativePermeability.scalers)
        props.RelativePermeability.scalers = ...
            getSubset(props.RelativePermeability.scalers, map.keep);
    end
    if isprop(props.RelativePermeability, 'swcon') ...
            && ~isempty(props.RelativePermeability.swcon)
        props.RelativePermeability.swcon = ...
            props.RelativePermeability.swcon(map.keep,:);
    end
end

%-------------------------------------------------------------------------%
function props = replacePVTProps(props, map)
    props = props.subset(map.keep);
end

%-------------------------------------------------------------------------%
function props = replaceFlowDisc(props, cellMap, faceMap, submodel)
    props = props.subset(cellMap.keep);
    if isprop(props.PhasePotentialDifference, 'pressureThreshold')
        props.PhasePotentialDifference.pressureThreshold = ...
            props.PhasePotentialDifference.pressureThreshold(faceMap.activeConn,:);
    end
    props = replaceOperators(props, submodel);
end


%-------------------------------------------------------------------------%
function discretization = replaceOperators(discretization, submodel)
    % Replace face upstream operator
    faceUpstr = submodel.operators.faceUpstr;
    fun       = @(sf) sf.UpwindDiscretization.setFunctionHandle(faceUpstr);
    names = {'FaceComponentMobility',   'FaceMobility'};
    discretization = replaceOperator(discretization, 'UpwindDiscretization', names, fun);
    % Replace kgrad operator (NB: Assuming simple two-point operator here!)
    tpfa = TwoPointFluxApproximation(submodel);
    fun  = @(sf) tpfa;
    names = {'PermeabilityPotentialGradient'};
    discretization = replaceOperator(discretization, 'PermeabilityGradientDiscretization', names, fun);
    % Replace gradient operator
    Grad = submodel.operators.Grad;
    fun  = @(sf) Grad;
    names = {'PressureGradient'};
    discretization = replaceOperator(discretization, 'Grad', names, fun);
end

%-------------------------------------------------------------------------%
function discretization = replaceOperator(discretization, operator, names, fun)
    names = discretization.getNamesOfStateFunctions()';
    for name = names
        sf = discretization.getStateFunction(name{1});
        if isprop(sf, operator)
            sf.(operator) = fun(sf);
            discretization = discretization.setStateFunction(name{1}, sf);
        end
    end
end

%-------------------------------------------------------------------------%
function prop = getSubset(prop, subs)
    fnames = reshape(fieldnames(prop), 1, []);
    for f = fnames
        if iscell(prop.(f{1}))
            prop.(f{1}) = cellfun(@(v) v(subs,:), prop.(f{1}), 'UniformOutput', false);
        elseif size(prop.(f{1}),1) == numel(subs)
            prop.(f{1}) = prop.(f{1})(subs,:);
        elseif isstruct(prop.(f{1}))
            prop.(f{1}) = getSubset(prop.(f{1}), subs);
        end
    end
end