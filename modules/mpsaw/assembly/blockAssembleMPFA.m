function assembly = blockAssembleMPFA(G, K, bcstruct, src, eta, globtbls, globmappings, varargin)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    opt = struct('verbose'  , mrstVerbose, ...
                 'blocksize', [], ...
                 'bcetazero', true, ...
                 'useVirtual', true);
    
    opt = merge_options(opt, varargin{:});
    
    blocksize = opt.blocksize;
    useVirtual = opt.useVirtual;
    
    nn = G.nodes.num;
    nblocks = floor(nn/blocksize);
    blocksizes = repmat(blocksize, nblocks, 1);
    if nn > nblocks*blocksize
        blocksizes = [blocksizes; nn - nblocks*blocksize];
    end
    nblocks = numel(blocksizes);
    blockinds = cumsum([1; blocksizes]);

    globcoltbl                = globtbls.coltbl;
    globcolrowtbl             = globtbls.colrowtbl;
    globcol2row2tbl           = globtbls.col2row2tbl;
    globnodetbl               = globtbls.nodetbl;
    globcellcoltbl            = globtbls.cellcoltbl;
    globnodecoltbl            = globtbls.nodecoltbl;
    globcellnodetbl           = globtbls.cellnodetbl;
    globcellnodecoltbl        = globtbls.cellnodecoltbl;
    globnodefacecoltbl        = globtbls.nodefacecoltbl;
    globcellcol2row2tbl       = globtbls.cellcol2row2tbl;
    globcellcolrowtbl         = globtbls.cellcolrowtbl;
    globcoltbl                = globtbls.coltbl;
    globcelltbl               = globtbls.celltbl;
    globnodetbl               = globtbls.nodetbl;
    globcellnodetbl           = globtbls.cellnodetbl;
    globnodefacetbl           = globtbls.nodefacetbl;
    globcellcoltbl            = globtbls.cellcoltbl;
    globnodecoltbl            = globtbls.nodecoltbl;
    globnodefacecoltbl        = globtbls.nodefacecoltbl;
    globcellnodefacetbl       = globtbls.cellnodefacetbl;
    globcellnodecoltbl        = globtbls.cellnodecoltbl;
    globcellnodecolrowtbl     = globtbls.cellnodecolrowtbl;
    globcellnodefacecolrowtbl = globtbls.cellnodefacecolrowtbl;
    globcolrowtbl             = globtbls.colrowtbl;
    globnodecolrowtbl         = globtbls.nodecolrowtbl;
    globcol2row2tbl           = globtbls.col2row2tbl;
    globcellcol2row2tbl       = globtbls.cellcol2row2tbl;
    globcellnodecol2row2tbl   = globtbls.cellnodecol2row2tbl;
    
    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = globnodefacetbl.num;
        globextflux = zeros(nf_num, 1);
    end

    if isempty(src)
        src = zeros(globcelltbl.num, 1);
    end
    
    globsrc = src;
    
    globbcdirichlet = bcstruct.bcdirichlet;
    globbcnodefacetbl = globbcdirichlet.bcnodefacetbl;
    globbcvals        = globbcdirichlet.bcvals;
    
    gnc = globcelltbl.num;
    gnbc = globbcnodefacetbl.num;
    
    globK = K;
    
    B11 = sparse(gnc, gnc);
    B12 = sparse(gnc, gnbc);
    B21 = sparse(gnbc, gnc);
    B22 = sparse(gnbc, gnbc);

    rhsc  = zeros(gnc, 1);
    rhsbc = zeros(gnbc, 1);
    
    for iblock = 1 : nblocks
        % Construction of tensor g (as defined in paper eq 4.1.2)
        nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';

        clear nodetbl;
        nodetbl.nodes = nodes;
        nodetbl = IndexArray(nodetbl);

        if opt.verbose
            fprintf('Assembling block %d/%d (%d nodes)\n', iblock, nblocks, nodetbl.num);
        end
        
        [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, 'useVirtual', useVirtual);

        celltbl     = tbls.celltbl;
        colrowtbl   = tbls.colrowtbl;
        nodetbl     = tbls.nodetbl;
        cellnodetbl = tbls.cellnodetbl;
        nodefacetbl = tbls.nodefacetbl;
        cellcoltbl  = tbls.cellcoltbl;
        nodecoltbl  = tbls.nodecoltbl;
        
        nodefacecoltbl        = tbls.nodefacecoltbl;
        cellnodefacetbl       = tbls.cellnodefacetbl;
        cellnodecoltbl        = tbls.cellnodecoltbl;
        cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
        cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
        cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
        nodecolrowtbl         = tbls.nodecolrowtbl;
        cellcol2row2tbl       = tbls.cellcol2row2tbl;
        cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;
        cellcolrowtbl         = tbls.cellcolrowtbl;
        
        globcell_from_cell = mappings.globcell_from_cell;
        
        map = TensorMap();
        map.fromTbl = globcellcolrowtbl;
        map.toTbl = cellcolrowtbl;
        map.mergefds = {'cells', 'coldim', 'rowdim'};
        
        map.pivottbl = cellcolrowtbl;
        c_num = celltbl.num; %shortcut
        ccr_num = cellcolrowtbl.num; %shortcut
        cr_num = colrowtbl.num; %shortcut
        gc_num = globcellcolrowtbl.num; %shortcut
        [cr, i] = ind2sub([cr_num, c_num], (1 : ccr_num)');
        map.dispind1 = sub2ind([cr_num, gc_num], cr, globcell_from_cell(i));
        map.dispind2 = (1 : ccr_num)';
        map.issetup = true;
        
        K = map.eval(globK);
        
        % We collect the degrees of freedom in the current block that belongs to the boundary.

        bcnodefacetbl = crossIndexArray(globbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        
        bcterm_exists = true;
        if bcnodefacetbl.num == 0
            bcterm_exists = false;
        end
        
        if bcterm_exists
            map = TensorMap();
            map.fromTbl = globbcnodefacetbl;
            map.toTbl = bcnodefacetbl;
            map.mergefds = {'nodes', 'faces'};
            bcind = map.getDispatchInd();
            
            clear bcdirichlet;
            bcdirichlet.bcnodefacetbl = bcnodefacetbl;
            bcdirichlet.bcvals = []; % not used locally
        else
            bcdirichlet = [];
        end

        % Setup main assembly matrices
        dooptimize = useVirtual;
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero, ...
                      'dooptimize', dooptimize);
        [matrices, bcvals, extra] = coreMpfaAssembly(G, K, bcdirichlet, tbls, mappings, opts);
    
        A11 = matrices.A11;
        A12 = matrices.A12;
        A21 = matrices.A21;
        A22 = matrices.A22;
        D = matrices.D;
        invA11 = matrices.invA11;
        
        map = TensorMap();
        map.fromTbl = globnodefacetbl;
        map.toTbl = nodefacetbl;
        map.mergefds = {'nodes', 'faces'};
        map = map.setup();
        
        extflux = map.eval(globextflux);

        locB11 = A22 - A21*invA11*A12;
        locrhsc = -A21*invA11*extflux; 
        
        if bcterm_exists
            locB12 = A21*invA11*D;
            locB21 = -D'*invA11*A12;
            locB22 = D'*invA11*D;
            locrhsbc = -D'*invA11*extflux;
        end

        % locB11 : celltbl       -> celltbl
        % locB22 : bcnodefacetbl -> bcnodefacetbl
        % locB12 : bcnodefacetbl -> celltbl
        % locB21 : celltbl       -> bcnodefacetbl 
        % We recall that these index arrays are all local.

        nc = celltbl.num;
        map = TensorMap();
        map.fromTbl = globcelltbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        cellind = map.getDispatchInd();
        
        [i, j, v] = find(locB11); 
        B11 = B11 + sparse(cellind(i), cellind(j), v, gnc, gnc);
        rhsc(cellind) = rhsc(cellind) + locrhsc;
        
        if bcterm_exists
            nbc = bcnodefacetbl.num;
            [i, j, v] = find(locB22);
            B22 = B22 + sparse(bcind(i), bcind(j), v, gnbc, gnbc);
            [i, j, v] = find(locB12);
            B12 = B12 + sparse(cellind(i), bcind(j), v, gnc, gnbc);
            [i, j, v] = find(locB21);
            B21 = B21 + sparse(bcind(i), cellind(j), v, gnbc, gnc);
            rhsbc(bcind) = rhsbc(bcind) + locrhsbc;
        end
        
    end
    
    B = [[B11, B12]; ...
         [B21, B22]];

    rhsc = rhsc + globsrc;
    rhsbc = rhsbc + globbcvals;
    
    rhs = [rhsc; ...
           rhsbc];
    
    assembly = struct( 'B'  , B  , ...
                       'rhs', rhs);    
end
