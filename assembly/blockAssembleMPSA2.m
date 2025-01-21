function assembly = blockAssembleMPSA2(G, prop, loadstruct, eta, globtbls, globmappings, varargin)
% Assembly of MPSA-weak
%
% Reference paper:
% Finite volume methods for elasticity with weak symmetry
% Keilegavlen, Eirik and Nordbotten, Jan Martin
% International Journal for Numerical Methods in Engineering
% 2017

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


    opt = struct('verbose'    , mrstVerbose, ...
                 'blocksize'  , []         , ...
                 'bcetazero'  , true       , ...
                 'useVirtual' , true       , ...
                 'extraoutput', false);
    opt = merge_options(opt, varargin{:});
    
    useVirtual = opt.useVirtual;
    blocksize  = opt.blocksize;
    
    nn = G.nodes.num;
    nblocks = floor(nn/blocksize);
    blocksizes = repmat(blocksize, nblocks, 1);
    if nn > nblocks*blocksize
        blocksizes = [blocksizes; nn - nblocks*blocksize];
    end
    nblocks = numel(blocksizes);
    blockinds = cumsum([1; blocksizes]);

    vectbl     = globtbls.vectbl;
    vec12tbl   = globtbls.vec12tbl;
    vec1212tbl = globtbls.vec1212tbl;
    
    globnodetbl        = globtbls.nodetbl;
    globfacetbl        = globtbls.facetbl;
    globnodefacetbl    = globtbls.nodefacetbl;
    globcellvectbl     = globtbls.cellvectbl;
    globnodevectbl     = globtbls.nodevectbl;
    globcellnodetbl    = globtbls.cellnodetbl;
    globcellnodevectbl = globtbls.cellnodevectbl;
    globnodefacevectbl = globtbls.nodefacevectbl;
    globcellvec1212tbl = globtbls.cellvec1212tbl;
    
    dim = vectbl.num;

    globextforce = loadstruct.extforce;
    globforce = loadstruct.force;

    globbc = loadstruct.bc;
    globbcnodefacetbl = globbc.bcnodefacetbl;        
    globbcnodefacetbl = globbcnodefacetbl.addLocInd('bcinds');    
    globbcnodefacevectbl = crossIndexArray(globbcnodefacetbl, vectbl, {}, ...
                                           'optpureproduct', true, 'virtual', useVirtual);
    globlinform = globbc.linform;
    globlinform = reshape(globlinform', [], 1);
    globlinformvals = globbc.linformvals;

    gncc = globcellvectbl.num;
    gnbc = globbcnodefacetbl.num;
    
    globC = setupStiffnessTensor2(prop, globtbls);
    
    B11 = sparse(gncc, gncc);
    B12 = sparse(gncc, gnbc);
    B21 = sparse(gnbc, gncc);
    B22 = sparse(gnbc, gnbc);

    rhscc = zeros(gncc, 1);
    rhsbc = zeros(gnbc, 1);
    
    % computes number of nodes per face
    
    map = TensorMap();
    map.fromTbl  = globnodefacetbl;
    map.toTbl    = globfacetbl;     
    map.mergefds = {'faces'};
    if useVirtual
        map.pivottbl = globnodefacetbl;
        map.dispind1 = (1 : globnodefacetbl.num)'
        map.dispind2 = globmappings.face_from_nodeface;
        map.issetup = true;
    else
        map = map.setup();
    end
    
    nnodesperface = map.eval(ones(globnodefacetbl.num, 1));

    for iblock = 1 : nblocks
        
        % Construction of tensor g (as defined in paper eq 4.1.2)
        nodes = [blockinds(iblock) : (blockinds(iblock + 1) - 1)]';

        clear nodetbl;
        nodetbl.nodes = nodes;
        nodetbl = IndexArray(nodetbl);

        if opt.verbose
            fprintf('Assembling block %d/%d (%d nodes)\n', iblock, nblocks, nodetbl.num);
        end
        
        [tbls, mappings] = setupStandardBlockTables2(G, nodetbl, globtbls, 'useVirtual', useVirtual, 'tblcase', 'mpsa');
        
        celltbl              = tbls.celltbl;
        facetbl              = tbls.facetbl;
        nodetbl              = tbls.nodetbl;
        cellnodetbl          = tbls.cellnodetbl;
        nodefacetbl          = tbls.nodefacetbl;
        cellvectbl           = tbls.cellvectbl;
        nodevectbl           = tbls.nodevectbl;
        nodefacevectbl       = tbls.nodefacevectbl;
        cellnodefacetbl      = tbls.cellnodefacetbl;
        cellnodevectbl       = tbls.cellnodevectbl;
        cellnodevec12tbl     = tbls.cellnodevec12tbl;
        cellnodefacevectbl   = tbls.cellnodefacevectbl;
        cellvec1212tbl       = tbls.cellvec1212tbl;

        globcell_from_cell         = mappings.globcell_from_cell;
        cell_from_cellnode         = mappings.cell_from_cellnode;
        node_from_cellnode         = mappings.node_from_cellnode;
        cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
        nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
        
        % Some shortcuts
        c_num     = celltbl.num;
        n_num     = nodetbl.num;
        cnf_num   = cellnodefacetbl.num;
        cnfc_num  = cellnodefacevectbl.num;
        cn_num    = cellnodetbl.num;
        cncr_num  = cellnodevec12tbl.num;
        nf_num    = nodefacetbl.num;
        nfc_num   = nodefacevectbl.num;
        d_num     = vectbl.num;
        
        dim = vectbl.num;
        
        % Obtain stiffness values for the block
        
        map = TensorMap();
        map.fromTbl  = globcellvec1212tbl;
        map.toTbl    = cellvec1212tbl;
        map.mergefds = {'cells', 'vec11', 'vec12', 'vec21', 'vec22'};
        
        map.pivottbl = cellvec1212tbl;
        cc2r2_num = cellvec1212tbl.num; %shortcut
        c2r2_num  = vec1212tbl.num; %shortcut
        gc_num    = globcellvec1212tbl.num; %shortcut
        [c2r2, i] = ind2sub([c2r2_num, c_num], (1 : cc2r2_num)');
        map.dispind1 = sub2ind([c2r2_num, gc_num], c2r2, globcell_from_cell(i));
        map.dispind2 = (1 : cc2r2_num)';
        map.issetup = true;
        
        C = map.eval(globC);
        
        % We collect the degrees of freedom in the current block that belongs to the boundary.
        
        [bcnodefacetbl, indstruct] = crossIndexArray(globbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        globbcnodeface_from_bcnodeface = indstruct{1}.inds;
        
        bcterm_exists = true;
        if bcnodefacetbl.num == 0
            bcterm_exists = false;
        end
        
        if bcterm_exists
            
            bcind = bcnodefacetbl.get('bcinds');
            bcnodefacevectbl = crossIndexArray(bcnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
            bcnodefacetbl = replacefield(bcnodefacetbl, {{'bcinds', ''}});
            
            linformvals = globlinformvals(bcind, :);

            map = TensorMap();
            map.fromTbl  = globbcnodefacevectbl;
            map.toTbl    = bcnodefacevectbl;
            map.mergefds = {'bcinds', 'vec', 'nodes', 'faces'};

            if useVirtual
                map.pivottbl = bcnodefacevectbl;
                [vec, i] = ind2sub([vectbl.num, bcnodefacetbl.num], (1 : bcnodefacevectbl.num)');
                map.dispind1 = sub2ind([vectbl.num, globbcnodefacetbl.num], vec, globbcnodeface_from_bcnodeface(i));
                map.dispind2 = (1 : bcnodefacevectbl.num)';
                map.issetup = true;
            else
                map = map.setup();
            end

            linform = map.eval(globlinform);
            linform = reshape(linform, dim, [])';
            
            bc = struct('bcnodefacetbl', bcnodefacetbl, ...
                        'linform'      , linform      , ...
                        'linformvals'  , linformvals);
        else
            
            bcnodefacetbl = replacefield(bcnodefacetbl, {{'bcinds', ''}});
            bc = struct('bcnodefacetbl', bcnodefacetbl, ...
                        'linform'      , []      , ...
                        'linformvals'  , []);
        end

        map = TensorMap();
        map.fromTbl  = globfacetbl;
        map.toTbl    = facetbl;
        map.mergefds = {'faces'};

        if useVirtual

            map.pivottbl = facetbl;
            map.dispind1 = mappings.globface_from_face;
            map.dispind2 = (1 : facetbl.num)';

            map.issetup = true;
            
        else
            
            map = map.setup();
            
        end

        nnpf = map.eval(nnodesperface);
        
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero);
        output = coreMpsaAssembly2(G, C, bc, nnpf, tbls, mappings, opts, 'useVirtual', useVirtual);

        matrices = output.matrices;
        bcvals   = output.bcvals;
        
        A11 = matrices.A11;
        A12 = matrices.A12;
        A21 = matrices.A21;
        A22 = matrices.A22;
        D = matrices.D;
        invA11 = matrices.invA11;

        % the solution is given by the system
        %
        % A = [[A11, A12, -D];
        %      [A21, A22,  0];
        %      [D' , 0  ,  0]];
        %
        % u = [u  (displacement at nodefacevectbl);
        %      u  (displacement at cellvectbl);
        %      lagmult];
        %
        % f = [extforce  (force at nodefacevectbl);
        %      force  (volumetric force at cellvectbl);
        %      bcvals (for the linear form at the boundary)];
        %
        % A*u = f
        %
        % Note: extforce is sparse and should only give contribution at facets
        % that are at the boundary
        %
        % By construction of the method, the matrix A11 is block-diagonal. Hence,
        % we invert it directly and reduce to a cell-centered scheme.


        % We reduced the system (shur complement) using invA11
        % We obtain system of the form
        %
        % B*u = rhs
        %
        % where
        %
        % B = [[B11, B12];
        %      [B21, B22]];
        %
        % u = [u (displacement at cellvectbl);
        %      lagmult];
        %
        % rhs = redextforce + [force;
        %                     bcvals]

        locB11 = A22 - A21*invA11*A12;
        
        if bcterm_exists
            locB12 = A21*invA11*D;
            locB21 = -D'*invA11*A12;
            locB22 = D'*invA11*D;
        end

        % locB11 : cellvectbl    -> cellvectbl
        % locB22 : bcnodefacetbl -> bcnodefacetbl
        % locB12 : bcnodefacetbl -> cellvectbl
        % locB21 : cellvectbl    -> bcnodefacetbl 
        % Above, we recall that these index arrays are all local.

        map = TensorMap();
        map.fromTbl  = cellvectbl;
        map.toTbl    = globcellvectbl;
        map.mergefds = {'cells', 'vec'};

        if useVirtual
            map.pivottbl = cellvectbl;

            [vec, i] = ind2sub([vectbl.num, celltbl.num], (1 : cellvectbl.num)');
            
            map.dispind1 = (1 : cellvectbl.num)';
            map.dispind2 = sub2ind([vectbl.num, tbls.globcelltbl.num], vec, globcell_from_cell(i));
            
            map.issetup = true;
            
        else
            map = map.setup();
        end
        % map.pivottbl will be cellvectbl so that dispind2 gives the index of
        % cellvectbl in globcellvectbl
        cellind = map.dispind2;
        
        ncc = cellvectbl.num;
        
        [i, j, v] = find(locB11);
        B11 = B11 + sparse(cellind(i), cellind(j), v, gncc, gncc);
        
        if bcterm_exists
            
            nbc = bcnodefacetbl.num;
            [i, j, v] = find(locB22);
            B22 = B22 + sparse(bcind(i), bcind(j), v, gnbc, gnbc);
            [i, j, v] = find(locB12);
            B12 = B12 + sparse(cellind(i), bcind(j), v, gncc, gnbc);
            [i, j, v] = find(locB21);
            B21 = B21 + sparse(bcind(i), cellind(j), v, gnbc, gncc);
            
        end
        
        % We get the part of external force that acts on the block
        map = TensorMap();
        map.fromTbl  = globnodefacevectbl;
        map.toTbl    = nodefacevectbl;
        map.mergefds = {'nodes', 'faces', 'vec'};

        if useVirtual
            
            map.pivottbl = nodefacevectbl;

            [vec, i] = ind2sub([vectbl.num, nodefacetbl.num], (1 : nodefacevectbl.num)');
            
            map.dispind1 = sub2ind([vectbl.num, globnodefacetbl.num], vec, mappings.globnodeface_from_nodeface(i));
            map.dispind2 = (1 : nodefacevectbl.num)';
            
            map.issetup = true;
        else
            map = map.setup();
        end

        extforce = map.eval(globextforce);
                
        % locrhscc in cellvectbl
        % locrhsbc in bcnodefacetbl

        locrhscc = -A21*invA11*extforce; 
        
        rhscc(cellind) = rhscc(cellind) + locrhscc;
        if bcterm_exists
            locrhsbc = -D'*invA11*extforce + bcvals;
            rhsbc(bcind) = rhsbc(bcind) + locrhsbc;
        end

    end
    
    B = [[B11, B12]; ...
         [B21, B22]];

    rhscc = rhscc + globforce;
    
    rhs = [rhscc; ...
           rhsbc];
    
    assembly = struct('B'       , B  , ...
                      'rhs'     , rhs, ...
                      'extforce', globextforce);

    if opt.extraoutput
        error('not yet implemented for block');
        assembly.divop = @(sol) mpsaDivOperator(sol, extforce, R1, R2, div);
    end
    
end
