function [substate0, substate, submodel, subforces, mappings] = buildSubProblem(model, state, state0, forces, subs)
%Undocumented Utility Function

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

    % Fint active cells and faces for subgrid----------------------
    G = model.G;
    Nc = G.cells.num;
    Nf = G.faces.num;

    op0 = model.operators;
    submodel = model;
    keepCells = false(Nc,1);
    keepCells(subs) = true;

    % Add all cells of wells intersecting the subset---------------
    for wNo = 1:numel(forces.W)
        wc = forces.W(wNo).cells;
        if any(keepCells(wc))
            keepCells(wc) = true;
        end
    end

    % Active cells
    active = keepCells(op0.N);
    activeConn = any(active,2);
    gc = op0.N(activeConn,:);
    gc = gc(~keepCells(gc));
    ghostCells = false(G.cells.num,1);
    ghostCells(gc) = true;
    keepCells(ghostCells) = true;

%             active = keepCells(op0.N);
%             activeConn = all(active,2);

    % Active faces
    cells = find(keepCells);
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1),1);
    keepFaces = false(G.faces.num,1);
    keepFaces(faces) = true;

    cells = find(keepCells & ~ghostCells);
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1),1);
    internalConn = false(G.faces.num,1);
    internalConn(faces) = true;
    bFaces = boundaryFaces(G);
    internalConn(bFaces) = false;
    internalConnParent = internalConn;
    internalConn = internalConn(keepFaces);            

    % Mappings between sub and full grid
    c_o2n = zeros(Nc, 1);
    c_o2n(keepCells) = 1:nnz(keepCells);
    c_n2o = find(keepCells);
    cellMap = struct('keep'   , keepCells  , ...
                     'old2new', c_o2n      , ...
                     'new2old', c_n2o      );

    f_o2n = zeros(Nf, 1);
    f_o2n(keepFaces) = 1:nnz(keepFaces);
    f_n2o = find(keepFaces);
    faceMap = struct('keep'   , keepFaces, ...
                     'old2new', f_o2n    , ...
                     'new2old', f_n2o    );

    mappings = struct('cellMap', cellMap, ...
                      'faceMap', faceMap);

    nf = nnz(activeConn);
    nc = nnz(keepCells);

    G.cells.ghost = false(G.cells.num,1);
    G.cells.ghost(ghostCells) = true;
    G.parent = G;

    G.cells.centroids = G.cells.centroids(keepCells, :);
    G.cells.volumes   = G.cells.volumes(keepCells);
    G.cells.diameters = G.cells.diameters(keepCells);
    G.cells.dx        = G.cells.dx(keepCells,:);
    G.cells.num = nc;
    G.cells.ghost = false(G.cells.num,1);
    G.cells.ghost(mappings.cellMap.old2new(ghostCells)) = true;

    G.faces.centroids = G.faces.centroids(keepFaces, :);
    G.faces.areas = G.faces.areas(keepFaces);
    if isfield(G.faces, 'dx')
        G.faces.diameters = G.faces.diameters(keepFaces,:);
        G.faces.dx = G.faces.dx(keepFaces,:);
    end
    G.faces.num = nnz(keepFaces);

    G.mappings = mappings;


    G.faces.neighbors = G.faces.neighbors(keepFaces,:);
    G.faces.neighbors(G.faces.neighbors~=0) = cellMap.old2new(G.faces.neighbors(G.faces.neighbors~=0));

    % Make operators for subgrid-----------------------------------
    T     = op0.T(activeConn);
    T_all = op0.T_all(keepFaces);
    N = cellMap.old2new(op0.N(activeConn, :));
    if size(N,2) == 1
        N = N';
    end
    M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
    C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);

    op = op0;
    op.C = C;
    op.M = M;
    upstr =  @(flag, x) faceUpstr(flag, x, N, [nf, nc]);
    op.faceUpstr = upstr;
    op.Div = @(x) C'*x;
    op.Grad = @(x) -C*x;
    op.faceAvg = @(x) M*x;
    op.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);

    disc = model.disc;
    vi = disc.velocityInterp;
    D = cell(1, G.griddim);
    for dNo = 1:G.griddim
        D{dNo} = vi.D{dNo}(keepCells, keepFaces);
    end
    vi.D = D;
    if G.griddim == 2
        vi.faceFlux2cellVelocity = @(v) [D{1}*v, D{2}*v];
    else
        vi.faceFlux2cellVelocity = @(v) [D{1}*v, D{2}*v, D{3}*v];
    end
    disc.velocityInterp = vi;

    op.T = T;
    op.T_all = T_all;
    op.N = N;
    op.pv = op0.pv(keepCells);
    op.internalConn = internalConn;

    % Copy state information---------------------------------------
    substate0 = state0;
    substate = state;
    flds = getCellFields();
    for fNo = 1:numel(flds)

        fn = flds{fNo};
        if isfield(substate, fn)

            ix = keepCells;
            if size(state.(fn)(:,1),1) > Nc
                ix = model.disc.getDofIx(state, Inf, cellMap.new2old);
            end
            substate.(fn) = state.(fn)(ix, :);

        end

        if isfield(substate0, fn)

            ix0 = keepCells;
            if size(state.(fn)(:,1),1) > Nc
                ix0 = model.disc.getDofIx(state0, Inf, cellMap.new2old);
            end
            substate0.(fn) = state0.(fn)(ix0, :);
        end

    end

    substate.flux = state.flux(keepFaces,:);

    % Make well struct for submodel--------------------------------
    W = forces.W;

    keepWells = false(numel(W), 1);
    keep = keepCells & ~ghostCells;
    for wNo = 1:numel(W)
%                 keepw = keepCells(W(wNo).cells);
        keepw = keep(W(wNo).cells);
        if any(keepw)
            assert(all(keepw));
            W(wNo).cells = cellMap.old2new(W(wNo).cells);
            keepWells(wNo) = true;
        end
    end
    W = W(keepWells);
    substate.wellSol  = substate.wellSol(keepWells);
    substate0.wellSol = substate0.wellSol(keepWells);

    % Make BC struct for submodel----------------------------------
    bc     = forces.bc;

    if ~isempty(bc)

        keepBC = any(bc.face == bFaces',2);
        bc.face       = faceMap.old2new(bc.face(keepBC));
        bc.type       = bc.type(keepBC);
        bc.value      = bc.value(keepBC);
        bc.sat      = bc.sat(keepBC,:);

    end

    % Gather forces in force struct--------------------------------
    subforces    = forces;
    subforces.W  = W;
    subforces.bc = bc;

    % Make subrock
    rock = submodel.rock;
    rock.poro = rock.poro(keepCells);
    rock.perm = rock.perm(keepCells,:);
    submodel.rock = rock;

    % Make submodel------------------------------------------------
    submodel.operators = op;
    submodel.G = G;
    submodel.FacilityModel = FacilityModel(submodel);
    submodel.FacilityModel = submodel.FacilityModel.setupWells(subforces.W);

    disc.G = G;
    disc.internalConn = op.internalConn;
    disc.internalConnParent = internalConnParent;
    disc.N = N;
    submodel.disc = disc;

    substate = submodel.disc.updateDofPos(substate);
    substate0 = submodel.disc.updateDofPos(substate0);

 end

 function flds = getCellFields()
    flds = {'pressure', 's', 'x', 'y', 'components', 'K', 'b', 'mob', ...
            'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w', 'bfactor', ...
            'w_p', 'dpressure', 'dpRel', 'dpAbs', 'switched', ...
            'nDof', 'degree', 'sdof'};
end
