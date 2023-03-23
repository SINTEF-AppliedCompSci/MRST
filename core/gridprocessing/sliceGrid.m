function [G, gix, G_slice] = sliceGrid(G, pnts, varargin)
% Generate cut-cell-grid (G) and lower dim grid G_slice from slices/cuts 
% defined by pnts and optional inputs 'normal' or 'cutDir' as described in 
% computeGridSlicePolygons.m. For multiple cuts/slices list (cell) of pnts
% and corresponding normal/cutDir are supported. See 
% computeGridSlicePolygons.m for further info on optional params.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('onlyFaces',  false, ...
             'topoSplit',  true, ...
             'mcompute',   exist('mcomputeGeometry', 'file') > 0, ...
             'gix_input',  []);
[opt, opt2] = merge_options(opt, varargin{:});
% others passed on to computeGridSlice

[pntsList, optList] = handleMultInput(pnts, opt2);
timer = tic();
npts  = numel(pntsList);

% establish cell-face-tags of new grid at the end, get map from original
% grid
faceTagMap = [];
if size(G.cells.faces, 2) == 2
    faceTagMap = getFaceTagMap(G);
end

for k = 1:npts
    % hard-set some opts that are required or not yet compatible with
    % grid-processing
    dispif(npts > 1 && mrstVerbose(), ...
        'Performing slice %d of %d\n', k, npts);

    poly = computeGridSlicePolygons(G, pntsList{k}, optList{k}{:}, ...
        'removePartialCuts'      , true, ...
        'outputForGridProcessing', true, ...
        'removeRepeated',          false);
    if G.griddim == 2
        [opt.topoSplit, opt.mcompute] = deal(false);
    end
    if k == 1
        gix = opt.gix_input;
    end
    if isempty(poly)
        warning('Slicing produced empty set, grid unchanged')
        G_slice = [];
        continue
    end
    if nargout >= 3 && G.griddim == 2
        warning('Output of 1D slice-grid not supported');
        G_slice = [];
    end
    
    [cutFaces, cutCells, coords, sliceFaces] = extractGridInfoFromPolygons(G, poly);
    
    [ncOld, nfOld] = deal(G.cells.num, G.faces.num);
    % don't update grid in case of no cutcells, but return gix, in case of
    % sliceFaces
    if ~isempty(cutCells.ix)
        % add coords to grid
        nNodesOrig = G.nodes.num;
        G.nodes.coords = [G.nodes.coords; coords];
        G.nodes.num    = size(G.nodes.coords, 1);
        % update node numbers. Deal with special case where node points to
        % grid-node (negative index)
        for kn = 1:size(cutFaces.nodes,2)
            nix = cutFaces.nodes(:,kn) < 0;
            cutFaces.nodes(nix, kn)  = -cutFaces.nodes(nix, kn);
            cutFaces.nodes(~nix, kn) =  cutFaces.nodes(~nix, kn) + nNodesOrig;
        end
        nix = cutCells.nodes < 0;
        cutCells.nodes(nix)  = -cutCells.nodes(nix);
        cutCells.nodes(~nix) =  cutCells.nodes(~nix) + nNodesOrig;
        
        % Perform actual grid splitting/cutting
        [G, isSplit, cutFaces, removedFaces]  = splitFaces(G, cutFaces, cutCells);
        if ~opt.onlyFaces
            [G, isCut, sliceFaces] = splitCells(G, cutCells, sliceFaces, opt.topoSplit);
        else
            isCut      = [];
            sliceFaces = [];
        end
        
        % reset type
        G.type = {'cutCellGrid'};
        % recompute geometry / bbox
        if opt.mcompute
            require libgeometry
            G = mcomputeGeometry(G);
        else
            G = computeGeometry(G);
        end
        G = addBoundingBoxFields(G);
        if any(G.cells.volumes == 0)
            nz = nnz(G.cells.volumes == 0);
            warning('Procedure produced %d cells with zero volume\n%s', nz, ...
                'Known issue: slicing precisely on boundary face');
        end
    else
        [isSplit, isCut] = deal([]);
    end
    gix = getIndices(ncOld, nfOld, cutCells.ix(isCut), cutFaces.ix(isSplit), sliceFaces, removedFaces, gix);
end
G = updateAdditionalFields(G, gix, faceTagMap);

dispif(mrstVerbose(), 'Performed %d slices in %fs\n', npts, toc(timer));
% create 2D slice-grid if requested
if nargout >= 3 && G.griddim == 3 && ~isempty(gix)
    G_slice = get2DGridFromFaces(G, gix.new.faces == 3);
end
end

%--------------------------------------------------------------------------
function [G, validIx, sliceFaces] = splitCells(G, cutCells, sliceFaces, topoSplit)
cix = cutCells.ix;
nc  = numel(cix);
% We store faces for the two new cells in cfo (replacing original) and cfn
%(new cell)
[cfo, cfn] = deal(cell(nc, 1)); % old/new cell-faces
nodePos    = cumsum([1;cutCells.nNodes]);
facePos    = G.cells.facePos;
validNum   = 1;
sliceIx    = logical(sparse(G.faces.num, 1));
failcount  = 0;
partcount  = 0;
topofail   = 0;
for k = 1:nc
    c = cix(k);
    % get original cell-faces
    fix = facePos(c):(facePos(c+1)-1);
    fo = G.cells.faces(fix);
    
    nix = nodePos(k):(nodePos(k+1)-1);
    ccn = cutCells.nodes(nix);
    
    if topoSplit && G.griddim == 3
        % topological splitting (more robust for curved faces)
        [np1, np2] = deal(G.faces.nodePos(fo), G.faces.nodePos(fo+1));
        fon  = G.faces.nodes(mcolon(np1, np2-1));
        [flag, success, isClean] = topologicalCellSplit(ccn, fon, cumsum([1;np2-np1]), G.faces.neighbors(fo)==c);
        % since new cell is given neighborposition 2, it has negative sign
        % according to node orientation
        flag = ~flag;
        if ~success
            topofail = topofail + 1;
        end
    else
        success = true;
    end
    if  ~topoSplit || ~success
        % split faces according to sign of
        % dot(plane normal, face-center - new-face-center)
        cent = mean(G.nodes.coords(ccn,:));
        v    = bsxfun(@minus, G.faces.centroids(fo,:), cent);
        %flag = v*normal < 0;
        flag = v*cutCells.normals(k,:)' < 0;
        success = true;
        isClean = true;
    end
    ok = (min(nnz(flag), nnz(~flag)) >= G.griddim) && success && isClean;
    
    if ok
        % order of cell-faces doesn't matter
        cfo{k} = [fo(flag); validNum + G.faces.num];
        cfn{k} = [fo(~flag); validNum + G.faces.num];
        
        % update neighbor-entries for neighboring faces of new cell
        nl = G.faces.neighbors(fo(~flag),:);
        nl(nl==c) = validNum + G.cells.num;
        G.faces.neighbors(fo(~flag),:) = nl;
        validNum = validNum +1;
    else
        % attempt to cut cell along face or partial cut
        % in case cut along face there should be a unique face with unique
        % flag
        if topoSplit && success
            if success
                if nnz(flag) >= G.griddim
                    flag = ~flag;
                end
                if nnz(flag) == 1
                    sliceIx(fo(flag)) = true;
                elseif ~isClean
                    % we have a partial cut, if this is to be allowed, the
                    % cut-face must be duplicated to get reasonnable topo
                    % and geometry.
                    partcount = partcount +1;
                else
                    failcount = failcount +1;
                end
            else
                failcount = failcount +1;
            end
        else
            ix = [];
            for kf = 1:numel(fo)
                npos = G.faces.nodePos(fo(kf) + [0;1]);
                nkf  = G.faces.nodes(npos(1):npos(2)-1);
                if numel(ccn) >= numel(nkf)
                    if all(ismember(nkf, ccn))
                        ix = kf;
                        break;
                    end
                end
            end
            if ~isempty(ix)
                sliceIx(fo(ix)) = true;
            else
                failcount = failcount +1;
            end
        end
    end
end
if any(sliceIx)
    sliceFaces = [sliceFaces; find(sliceIx)];
end
if failcount > 0
    warning('Cutting of %d out of %d cells failed', failcount, nc);
end
if partcount >0
    dispif(mrstVerbose, 'Ignored processing of %d partially cut cells \n', partcount);
end
if topofail >0
    dispif(mrstVerbose, 'Topological split failed for %d out of %d cells \n', topofail, nc);
end
% include info for new faces
validIx = find(~cellfun(@isempty, cfn));
nc = numel(validIx);
G.faces.neighbors = [G.faces.neighbors; [cix(validIx), (1:nc)'+G.cells.num]];
G.faces.num       = G.faces.num + nc;
G.faces.tag       = [G.faces.tag; zeros(nc,1)];
if numel(validIx) < numel(cfn)
    nnpos = cumsum([1;cutCells.nNodes]);
    ii    = mcolon(nnpos(validIx), nnpos(validIx+1)-1);
    G.faces.nodes     = [G.faces.nodes; cutCells.nodes(ii)];
else
    G.faces.nodes     = [G.faces.nodes; cutCells.nodes];
end
G.faces.nodePos   = [G.faces.nodePos; G.faces.nodePos(end) + cumsum(cutCells.nNodes(validIx))];

G.cells.num = G.cells.num + nc;
% use special purpose function to update facePos/faces
[f1, fp1] = replaceEntries(G.cells.faces, cfo(validIx), cix(validIx), G.cells.facePos);
G.cells.facePos = [fp1; fp1(end) + cumsum(cellfun(@numel, cfn(validIx)))];
G.cells.faces   = [f1; vertcat(cfn{validIx})];
end

%--------------------------------------------------------------------------
function [G, validIx, cutFaces, removedFaces] = splitFaces(G, cutFaces, cutCells)
[fix, cix] = deal(cutFaces.ix, cutCells.ix);
[nf , nc ] = deal(numel(fix), numel(cix));
% Need to update original and create new face-nodes for cutfaces
[fno, fnn] = deal(cell(nf,1));
% Need to update original cell-faces for cutcells
cfo   = cell(nc, 1);
cfpos = G.cells.facePos;
for kc = 1:numel(cfo)
    ic = cfpos(cix(kc)):(cfpos(cix(kc)+1)-1);
    cfo{kc} = G.cells.faces(ic, 1);
end
% Need to account for that additonal (non-cutcell) neighboring cells are
% influenced. At most one per face.
[cixe, cfe] = deal(cell(nf,1));
% create look-up index for cut-cells
cutCellIx  = zeros(G.cells.num,1);
cutCellIx(cix) = (1:nc)';
% centroids for new faces
centN = nan(nf,3);
nodePos = G.faces.nodePos;
validCount = 0;
for k = 1:nf
    f = fix(k);
    % get original face-nodes
    nix = nodePos(f):(nodePos(f+1)-1);
    no  = G.faces.nodes(nix);
    nn  = cutFaces.nodes(k,:);
    p   = cutFaces.nodePos(k,:);
    if G.griddim == 3
        if p(1) > p(2)
            p  = p([2,1]);
            nn = nn([2,1]);
        end
        % redistribute nodes according to given position (p)
        ix1 = 1:p(1);
        ix2 = (p(1)+1):p(2);
        ix3 = (p(2)+1):numel(no);
        if cutFaces.isTrueSplit(k)
            fno{k} = [nn(1); no(ix2); nn(2)];
            fnn{k} = [nn(2); no(ix3); no(ix1); nn(1)];
        else % just add nodes to current face
            fno{k} = [no(ix1); nn(1); no(ix2); nn(2); no(ix3)];
        end
    elseif G.griddim == 2
        if cutFaces.isTrueSplit(k)
            fno{k} = [no(1); nn];
            fnn{k} = [nn;    no(2)];
        else
            fno{k} = [no(1); nn; no(2)];
        end
    end
    % if new nodes hit grid-nodes, we get repeated fno, fix this
    validNewFace = true;
    isSame = [false; diff(fno{k}) == 0];
    if any(isSame)
        fno{k}(isSame) = [];
    end
    isSame = [false; diff(fnn{k}) == 0];
    if any(isSame)
        fnn{k}(isSame) = [];
    end
    if numel(fno{k}) < G.griddim
        [fno{k}, fnn{k}] = deal(fnn{k}, fno{k});
    end
    if numel(fnn{k}) < G.griddim
        fnn{k} = [];
        validNewFace = false;
    end
    validCount = validCount + validNewFace;
    % temporary update face-center (used in splitCells)
    G.faces.centroids(f,:) = mean(G.nodes.coords(fno{k},:));
    if validNewFace
        centN(k,:) = mean(G.nodes.coords(fnn{k},:));
        % neighboring cells need to update their cellfaces
        c = G.faces.neighbors(f,:);
        c = c(c>0);
        for kc = 1:numel(c)
            ccix = cutCellIx(c(kc));
            if ccix > 0 % cut-cell
                cfo{ccix} = [cfo{ccix}; validCount+G.faces.num];
            else % we have a cell that is not a cut-cell
                ic = cfpos(c(kc)):(cfpos(c(kc)+1)-1);
                cfe{k}  = [G.cells.faces(ic,1); validCount+G.faces.num];
                cixe{k} = c(kc);
            end
        end
    end
end
% reduce list of influenced 'non-cutcell' cells/cellfaces
ri = cellfun(@isempty, cixe);
[cfe, cixe] = deal(cfe(~ri), vertcat(cixe{~ri}));
% check for repeated cell-indices, could happen for partly collapsed cells
% identify candidate pairs of faces that might be duplicate/collapsed
[~, ia, ib] = unique(cixe, 'stable');
duplicatesExist = false;
if numel(ia) ~= numel(ib)
    duplicatesExist = true;
    [~, occurance] = rlencode(ib);
    if any(occurance>2)
        warning('Not able to resolve potential duplicate faces');
    end
    ii = find(occurance==2);
    duplicateCandidates = repmat({nan(numel(ii), 2)}, [1,2]);
    adjecentFaces = fix(~ri);
    for k = 1:numel(ii)
        cc = find(ib==ii(k));
        % include all potential faces, remove duplicate later
        cfNew = union(cfe{cc(1)}, cfe{cc(2)});
        duplicateCandidates{1}(k,:) = reshape(adjecentFaces(cc), [1,2]);
        duplicateCandidates{2}(k,:) = setxor(cfe{cc(1)}, cfe{cc(2)});
        % set both, first occuring will be overwritten
        [cfe{cc}] = deal(cfNew);
    end
end
    
validIx = find(~cellfun(@isempty, fnn));
nf = numel(validIx);

G.faces.num       = G.faces.num + nf;
% neighbors of new faces is just a copy of corresponding for original
G.faces.neighbors = [G.faces.neighbors; G.faces.neighbors(fix(validIx),:)];
G.faces.tag       = [G.faces.tag; zeros(nf,1)];
% temporary centroids
G.faces.centroids = [G.faces.centroids; centN(validIx,:)];
% update facenodes/nodepos
%[n1, np1] = replaceEntries(G.faces.nodes, fno(validIx), fix(validIx), G.faces.nodePos);
[n1, np1] = replaceEntries(G.faces.nodes, fno, fix, G.faces.nodePos);
G.faces.nodePos = [np1; np1(end) + cumsum(cellfun(@numel, fnn(validIx)))];
G.faces.nodes   = [n1; vertcat(fnn{validIx})];
% update cellfaces/facepos for cutcells and potenial extras
[G.cells.faces, G.cells.facePos] = ...
    replaceEntries(G.cells.faces(:,1), [cfo;cfe], [cix; cixe], G.cells.facePos);
removedFaces = [];
if duplicatesExist
   [G, removedFaces] = removeDuplicateFaces(G, duplicateCandidates);
end
end

%--------------------------------------------------------------------------
function [flag, success, isClean] = topologicalCellSplit(n1, n2, pos, isOut)
% comparing node-lists for two adjacent faces, there should be at least two
% common (potentially more as a result of snapping ?) adjecant nodes, set
% flag according to
% true  : node order same and isOut/node order reverse and ~isOut
% if face is not ajecent to n1 it gets the same sign as any other of its
% ajecent faces (iteratively)
success = true;
nAdj  = zeros(size(n1));
nf      = numel(pos)-1;
flag    = false(nf, 1);
found   = false(nf, 1);
for k = 1:nf
    nb = n2(pos(k):pos(k+1)-1);
    [~, ia, ib] = intersect(n1, nb, 'stable');
    if numel(ia) >= 2
        if numel(ia) == 2
            startIx = 2 - (diff(ia)==1);
            chkIx = 1:2;
        else
            % special case might occur for overlapping faces
            startIx = [abs(diff(ia))==1; ia(1)==1 && ia(end)==numel(n1)];
            chkIx = find(startIx, 1, 'first') + (0:1);    
            %ii = find(abs(diff(ia))~=1);
            %if isempty(ii)
            %    startIx = 1:numel(ia)-1;
            %else
            %    startIx = mod(ii+(1:numel(ia)-1), ;
            %end
        end
        nAdj(ia(startIx)) = nAdj(ia(startIx)) + 1;
        s1 = mod(diff(ia(chkIx)), numel(n1)) == 1;
        s2 = mod(diff(ib(chkIx)), numel(nb)) == 1;
        flag(k)  = (s1==s2) == isOut(k);
        found(k) = true;
    end
end
% for a clean cut, each cut-cell segment must have two adjecent faces
isClean = all(nAdj == 2);
if any(~found) % there are faces not adjecent to cut-cell face
    its = 0;
    maxIts = 100; % should be plenty
    while any(~found) && its < maxIts
        ix = find(~found);
        for k = 1:numel(ix)
            nix = n2(pos(ix(k)):pos(ix(k)+1)-1);
            for fk = 1:numel(flag)
                if found(fk)
                    nk = n2(pos(fk):pos(fk+1)-1);
                    if nnz(ismember(nix, nk)) >= 2
                        flag(ix(k))  = flag(fk);
                        found(ix(k)) = true;
                        break;
                    end
                end
            end
        end
        its = its +1;
    end
    if its == maxIts
        success = false;
    end
end
end

%--------------------------------------------------------------------------
function [G, removedFaces] = removeDuplicateFaces(G, fcand)
% Partly collpased cells that are cut might result in a degenerate cut
% polygon, this function identifies pairs of collapsed faces from a list 
% of candidates and removes them from the cell that contains both
nd = size(fcand{1},1);
isDuplicate = true(nd,1);
fdup = nan(nd, 2);
% find duplicates from candidates
for k = 1:nd
    for kc = 1:2
        [f1, f2] = deal(fcand{kc}(k,1), fcand{kc}(k,2));
        np1  = G.faces.nodePos([f1, f1+1]) - [0 1]';
        nix1 = G.faces.nodes(np1(1):np1(2));
        np2 = G.faces.nodePos([f2, f2+1]) - [0 1]';
        nix2 = G.faces.nodes(np2(1):np2(2));
        if ~all(sort(nix1) == sort(nix2) )
            isDuplicate(k) = false;
%             Might want to check coordinates for non-unique nodes?
%             cc1 = sortrows(G.nodes.coords(nix1,:));
%             cc2 = sortrows(G.nodes.coords(nix2,:));
%             if ~all(all(abs(cc1-cc2)<sqrt(eps)))
%                 isDuplicate(k) = false;
%             end
        else
            isDuplicate(k) = true;
            fdup(k,:) = [f1, f2];
            break;
        end
    end
end
if ~all(isDuplicate)
    warning('Potential problem with handling duplicate faces')
end
f = fdup(isDuplicate,:);
nd = size(f,1);
cfNew = cell(nd,1);
cellUpdate = nan(nd, 1);
for k = 1:nd
    [f1, f2] = deal(f(k,1), f(k,2));
    n1 = G.faces.neighbors(f1,:);
    n2 = G.faces.neighbors(f2,:);
    common = intersect(n1, n2);
    assert(numel(common)==1);
    n1(n1 == common) = n2(n2 ~= common);
    G.faces.neighbors(f1, :) = n1;
    G.faces.neighbors(f2, :) = [nan, nan];
    % remove both collapsed faces from cell
    fpos = G.cells.facePos(common:(common+1)) - [0 1]';
    cfNew{k} = G.cells.faces(fpos(1):fpos(2));
    cfNew{k} = setdiff(cfNew{k}, [f1; f2]);
    cellUpdate(k) = common;
end
keepIx = true(G.faces.num,1);
keepIx(f(:,2)) = false;
reindex = cumsum(keepIx);
reindex(f(:,2)) = reindex(f(:,1));
fldnm = fieldnames(G.faces);
for k = 1:numel(fldnm)
    nm = fldnm{k};
    if size(G.faces.(nm), 1) == G.faces.num
        G.faces.(nm) = G.faces.(nm)(keepIx,:);
    end
end
G.faces.num = G.faces.num - nd;
[cellFaces, G.cells.facePos] = ...
    replaceEntries(G.cells.faces, cfNew, cellUpdate, G.cells.facePos);
G.cells.faces = reindex(cellFaces);

emptyFaceNodes = cell(size(f,1), 1);
[G.faces.nodes, G.faces.nodePos] = ...
    replaceEntries(G.faces.nodes, emptyFaceNodes, f(:,2), G.faces.nodePos);
removedFaces = f(:,2);
end

%--------------------------------------------------------------------------
function ix = getIndices(nc, nf, cutCellIx, cutFaceIx, sliceFaces, removedFaces, ixp)
facesPrev = (1:nf)'; 
if any(removedFaces)
    removedPrev = removedFaces(removedFaces <= nf);
    facesPrev(removedPrev) = [];
    removedCur = removedFaces(removedFaces > nf);
    cutFaceIx(removedCur - nf) = [];
    nf = numel(facesPrev);
end
[nc2, nf2] = deal(numel(cutCellIx), numel(cutFaceIx));
if isempty(ixp)
    [nc_orig, nf_orig] = deal(nc, nf);
    ix.old    = struct('cells', true(nc,1), 'faces', true(nf,1));
    ix.new    = struct('cells', ones(nc+nc2,1), 'faces', ones(nf+nf2+nc2,1));
    ix.parent = struct('cells', [(1:nc)'; cutCellIx], 'faces', [facesPrev; cutFaceIx; zeros(nc2,1)]);
    previous  = ix.parent;
else
    if any(removedFaces)
        ixp.parent.faces(facesPrev) = [];
        ixp.new.faces(facesPrev) = [];
    end
    [nc_orig, nf_orig] = deal(numel(ixp.old.cells), numel(ixp.old.faces));
    ix.old    = ixp.old;
    ix.new    = struct('cells', [ixp.new.cells; ones(nc2, 1)], ...
                       'faces', [ixp.new.faces; ones(nf2+nc2, 1)]);
    ix.parent = struct('cells', [ixp.parent.cells; ixp.parent.cells(cutCellIx)], ...
                       'faces', [ixp.parent.faces; ixp.parent.faces(cutFaceIx); zeros(nc2,1)]);
    previous = struct('cells', [(1:nc)'; cutCellIx], 'faces', [facesPrev; cutFaceIx; zeros(nc2,1)]);
end
% set split cells/faces to false
ix.old.cells(cutCellIx(cutCellIx < nc_orig)) = false;
ix.old.faces(cutFaceIx(cutFaceIx < nf_orig)) = false;
% set new/updated cells to 2
ix.new.cells([cutCellIx; nc + (1:nc2)']) = 2;
% set new/updated faces to 2, faces along plane to 3
ix.new.faces([cutFaceIx; nf + (1:nf2)']) = 2;
ix.new.faces(nf + nf2 + (1:nc2)) = 3;
ix.new.faces(sliceFaces) = 3;
if ~isempty(ixp)
    %pfac   = ix.parent.faces;
    pfac = previous.faces;
    pfacix = pfac > 0; 
    ix.new.faces(pfacix) = max(ix.new.faces(pfacix), ixp.new.faces(pfac(pfacix)));
    %pcel   = ix.parent.cells;
    pcel   = previous.cells;
    pcelix = pcel > 0;
    ix.new.cells(pcelix) = max(ix.new.cells(pcelix), ixp.new.cells(pcel(pcelix)));
end
end

%--------------------------------------------------------------------------
function [w, posNew] = replaceEntries(v, vi, ix, pos)
% helper for replacing cell-face/face-node entries
ni = cellfun(@numel, vi);
nNew = deal(diff(pos));
nNew(ix) = ni;
posNew = cumsum([1; nNew]);
nw = posNew(end)-1;
nb = numel(posNew)-1;
% make linear index to ~ix
nix = true(nb, 1);
nix(ix) = false;
nix = find(nix);
ix1 = mcolon(pos(nix), pos(nix+1)-1);
ix2 = mcolon(posNew(nix), posNew(nix+1)-1);
ixn = mcolon(posNew(ix),  posNew(ix+1)-1);
w  = zeros(nw,1);
w(ix2) = v(ix1);
w(ixn) = vertcat(vi{:});
% if set is empty, reduce posNew
remove = cellfun(@isempty, vi);
posNew(ix(remove)) = [];
end

%--------------------------------------------------------------------------
function faceTagMap = getFaceTagMap(G)
nc  = G.cells.num;
cellNo = rldecode((1:nc)', diff(G.cells.facePos));
isPos  = G.faces.neighbors(G.cells.faces(:,1),1) == cellNo;
cfTag  = G.cells.faces(:,2);
[faceTagMap.pos, faceTagMap.neg] = deal(zeros(G.faces.num, 1));
faceTagMap.pos(G.cells.faces(isPos))  = cfTag(isPos);
faceTagMap.neg(G.cells.faces(~isPos)) = cfTag(~isPos);
end

%--------------------------------------------------------------------------
function G = updateAdditionalFields(G, gix, faceTagMap)
if ~isempty(faceTagMap)
    nc  = G.cells.num;
    ncf = size(G.cells.faces,1);
    cellNo = rldecode((1:nc)', diff(G.cells.facePos));
    isPos  = G.faces.neighbors(G.cells.faces(:,1),1) == cellNo;
    cfTag   = zeros(ncf,1);
    fParent = gix.parent.faces(G.cells.faces(:,1));
    nzero   = fParent > 0;
    cfTag(isPos & nzero)  = faceTagMap.pos(fParent(isPos & nzero));
    cfTag(~isPos & nzero) = faceTagMap.neg(fParent(~isPos & nzero));
    G.cells.faces = [G.cells.faces(:,1), cfTag];
end
if isfield(G.cells, 'indexMap') && numel(G.cells.indexMap) < G.cells.num
    G.cells.indexMap = G.cells.indexMap(gix.parent.cells);
end
end

%--------------------------------------------------------------------------
function [p, opt] = handleMultInput(p, opt)
hasNormal = true;
dpos = find(strcmp(opt(1:2:end), 'cutDir'));
if ~isempty(dpos)
    v = opt{2*dpos};
    hasNormal = false;
end
npos = find(strcmp(opt(1:2:end), 'normal'));
if ~isempty(npos)
    if ~hasNormal
        warning('Ignoring optionial input ''normal'' since non-empty cut-direction');
    else
        v = opt{2*npos};
        hasNormal = true;
    end
end
hasRadii = false;
rpos     = find(strcmp(opt(1:2:end), 'radius'));
if ~isempty(rpos)
    r = opt{2*rpos};
    hasRadii = size(r,1) > 1;
    if hasRadii && ~iscell(r)
       %r = mat2cell(r, [ones(numel(r),size(r,2)), 1]);
       r = mat2cell(r, ones(numel(r),size(r,2)), 1);
    end
end 
if ~iscell(v)
    nv = size(v,1);
    v  = mat2cell(v, ones(nv,1), 3);
end
if ~iscell(p)
    np = size(p,1);
    if hasNormal % individual point
        p = mat2cell(p, ones(np,1), 3);
    else
        p = {p};
    end
end
[np, nv] = deal(numel(p), numel(v));
num = max([np, nv]);
if num > 1
    if np == 1
        p = repmat(p, [num, 1]);
    elseif nv == 1
        v = repmat(v, [num, 1]);
    end
end
[np, nv] = deal(numel(p), numel(v));
assert(np==nv && (~hasRadii || np == numel(r)), 'Can''t work out input');
optlist = repmat({opt}, [num, 1]);
for k = 1:num
    if ~hasNormal
        pos = dpos;
    else
        pos = npos;
    end
    optlist{k}{2*pos} = v{k};
    if hasRadii
        optlist{k}{2*rpos} = r{k};
    end
end
opt = optlist;
end
