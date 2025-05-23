function assembly = blockAssembleBiot(G, props, drivingforces, eta, globtbls, globmappings, varargin)
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


    opt = struct('verbose'         , mrstVerbose, ...
                 'assemblyMatrices', false, ...
                 'addAdOperators'  , false, ...
                 'blocksize'       , []   , ...
                 'bcetazero'       , true , ...
                 'useVirtual'      , false , ...
                 'extraoutput'     , false);
        
    opt = merge_options(opt, varargin{:});

    % option 'addAdOperators' not supported for the moment
    assert(~opt.addAdOperators, 'option not implemented yet');
    
    % We solve the system
    %
    %  A*u = f
    %
    % where
    %
    %
    %        | A11    A12    0      A14    A15    0    |
    %        | A21    A22    0      0      0      0    |
    %  A =   | 0      0      A33    A34    0      A36  |
    %        | A41    A42    A43    A44    0      0    |
    %        | A51    0      0      0      0      0    |
    %        | 0      0      A63    0      0      0    | 
    %
    %
    %       | displacement_nfc (node face dofs, belongs to nodefacevectbl)         |
    %       | displacement_c   (cell dofs, belongs to cellvectbl)                  |
    %  u =  | pressure_nf      (node face dofs, belongs to nodefacetbl)            |
    %       | pressure_c       (cell dofs, belongs to celltbl)                     |
    %       | lambda1          (lagrangian multiplier for Dirichlet mechanical bc) |
    %       | lambda2          (lagrangian multiplier for Dirichlet fluid bc) |
    %
    %
    %       | exterior forces      |
    %       | volumetric forces    |
    %  f =  | exterior fluxes      |
    %       | source               |
    %       | mechanical bc values |              
    %       | fluid bc values      |

    
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

    vectbl             = globtbls.vectbl;
    vec12tbl           = globtbls.vec12tbl;
    vec1212tbl         = globtbls.vec1212tbl;
    globcellvectbl     = globtbls.cellvectbl;
    globcellnodetbl    = globtbls.cellnodetbl;
    globnodefacevectbl = globtbls.nodefacevectbl;
    globcellvec12tbl   = globtbls.cellvec12tbl;
    globcelltbl        = globtbls.celltbl;
    globfacetbl        = globtbls.facetbl;
    globnodefacetbl    = globtbls.nodefacetbl;
    globcellvec1212tbl = globtbls.cellvec1212tbl;
    
    dim = vectbl.num;
    
    mechprops  = props.mechprops;
    fluidprops = props.fluidprops;
    coupprops  = props.coupprops;
   
    globC     = setupStiffnessTensor(mechprops, globtbls, globmappings, 'useVirtual', useVirtual);
    globK     = fluidprops.K;
    globalpha = coupprops.alpha;
    globrho   = coupprops.rho;
    
    loadstruct = drivingforces.mechanics;
    fluidforces = drivingforces.fluid;
    bcstruct = fluidforces.bcstruct;
    src = fluidforces.src;
    
    % setup global structure for mechanic boundary condition
    globextforce = loadstruct.extforce;
    globforce = loadstruct.force;

    globmechbc = loadstruct.bc;
    globmechbcnodefacetbl = globmechbc.bcnodefacetbl;        
    globmechbcnodefacetbl = globmechbcnodefacetbl.addLocInd('bcinds');    
    globmechbcnodefacevectbl = crossIndexArray(globmechbcnodefacetbl, vectbl, {}, ...
                                           'optpureproduct', true);
    globlinform = globmechbc.linform;
    globlinform = reshape(globlinform', [], 1);
    globlinformvals = globmechbc.linformvals;
    
    % setup global structure for fluid boundary condition (globextflux and globsrc)
    if ~isempty(bcstruct.bcneumann)
        error('not yet implemented');
    else
        nf_num = globnodefacetbl.num;
        globextflux = zeros(nf_num, 1);
    end

    if isempty(src)
        globsrc = zeros(globcelltbl.num, 1);
    else
        globsrc = src;
    end
    
    globbcdirichlet = bcstruct.bcdirichlet;
    globfluidbcnodefacetbl = globbcdirichlet.bcnodefacetbl;
    globfluidbcvals = globbcdirichlet.bcvals;

    % setup assembly matrices where the values computed at block level will be stored
    %       | B11  B12  B13      |
    %  B =  | B21  B22  B23  B24 |
    %       | B31  B32  B33      |
    %       |      B42       B44 |
    %
    %
    %       | displacement at cell          |    in globcellvectbl
    %  u =  | pressure at cell              |    in globcelltbl
    %       | lagrange multiplier mechanics |    in globmechbcnodefacetbl
    %       | lagrange multiplier fluid     |    in globfluidbcnodefacetbl
    

    % Setup the global matrices B  
    bglobtbls = {globcellvectbl, globcelltbl, globmechbcnodefacetbl, globfluidbcnodefacetbl};
    nA = 6;
    nB = 4;
    B = cell(nB, 1);
    rhs = cell(nB, 1);
    bglobnums = cell(nB, 1);
    
    for i = 1 : nB
        B{i} = cell(nB, 1);
        insB{i} = cell(nB, 1);
        ni = bglobtbls{i}.num;
        bglobnums{i} = ni;
        rhs{i} = zeros(ni, 1);
        for j = 1 : nB;
            nj = bglobtbls{j}.num;
            B{i}{j} = sparse(ni, nj);
        end
    end
    
    if opt.assemblyMatrices
        % for debugging: store the non-reduced system
        aglobtbls = {globnodefacevectbl, globcellvectbl, globnodefacetbl, globcelltbl, globmechbcnodefacetbl, globfluidbcnodefacetbl};
        A = cell(nA, 1);
        for i = 1 : nA
            A{i} = cell(nA, 1);
            insA{i} = cell(nA, 1);
            ni = aglobtbls{i}.num;
            for j = 1 : nA
                nj = aglobtbls{j}.num;
                A{i}{j} = sparse(ni, nj);
                insA{i}{j} = sparse(ni, nj);
            end
        end
    end

    % we compute the number of nodes per cells and per faces
    map = TensorMap();
    map.fromTbl  = globcellnodetbl;
    map.toTbl    = globcelltbl;     
    map.mergefds = {'cells'};
    map = map.setup();
    
    nnodespercell = map.eval(ones(globcellnodetbl.num, 1));
    
    map = TensorMap();
    map.fromTbl  = globnodefacetbl;
    map.toTbl    = globfacetbl;     
    map.mergefds = {'faces'};
    map = map.setup();
    
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
        
        [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, 'useVirtual', useVirtual, 'tblcase', 'both');

        celltbl     = tbls.celltbl;
        vec12tbl    = tbls.vec12tbl;
        nodetbl     = tbls.nodetbl;
        facetbl     = tbls.facetbl;
        cellnodetbl = tbls.cellnodetbl;
        nodefacetbl = tbls.nodefacetbl;
        cellvectbl  = tbls.cellvectbl;
        nodevectbl  = tbls.nodevectbl;
        
        nodefacevectbl     = tbls.nodefacevectbl;
        cellnodefacetbl    = tbls.cellnodefacetbl;
        cellnodevectbl     = tbls.cellnodevectbl;
        cellnodevec12tbl   = tbls.cellnodevec12tbl;
        cellnodefacevectbl = tbls.cellnodefacevectbl;
        cellvec1212tbl     = tbls.cellvec1212tbl;
        cellvec12tbl       = tbls.cellvec12tbl;
        
        globcell_from_cell = mappings.globcell_from_cell;
    
        % Obtain stiffness values for the block
        map = TensorMap();
        map.fromTbl  = globcellvec1212tbl;
        map.toTbl    = cellvec1212tbl;
        map.mergefds = {'cells', 'vec11', 'vec12', 'vec21', 'vec22'};
        
        map.pivottbl = cellvec1212tbl;
        c_num     = celltbl.num;             % shortcut
        gc_num    = globcellvec1212tbl.num; % shortcut
        cc2r2_num = cellvec1212tbl.num;     % shortcut
        c2r2_num  = vec1212tbl.num;         % shortcut
        [c2r2, i] = ind2sub([c2r2_num, c_num], (1 : cc2r2_num)');
        map.dispind1 = sub2ind([c2r2_num, gc_num], c2r2, globcell_from_cell(i));
        map.dispind2 = (1 : cc2r2_num)';
        map.issetup = true;
        
        C = map.eval(globC);

        % Obtain permeability values for the block
        map = TensorMap();
        map.fromTbl  = globcellvec12tbl;
        map.toTbl    = cellvec12tbl;
        map.mergefds = {'cells', 'vec1', 'vec2'};
        
        map.pivottbl = cellvec12tbl;
        ccr_num = cellvec12tbl.num; % shortcut
        cr_num = vec12tbl.num;      % shortcut
        [cr, i] = ind2sub([cr_num, c_num], (1 : ccr_num)');
        map.dispind1 = sub2ind([cr_num, gc_num], cr, globcell_from_cell(i));
        map.dispind2 = (1 : ccr_num)';
        map.issetup = true;
        
        K = map.eval(globK);

        
        % Assembly mechanical part
        
        % We collect the degrees of freedom in the current block that belongs to the boundary.
        [mechbcnodefacetbl, indstruct] = crossIndexArray(globmechbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        globmechbcnodeface_from_mechbcnodeface = indstruct{1}.inds;
        
        mechbcterm_exists = true;
        if mechbcnodefacetbl.num == 0
            mechbcterm_exists = false;
        end
                
        if mechbcterm_exists
            
            mechbcinds = mechbcnodefacetbl.get('bcinds');
            mechbcnodefacevectbl = crossIndexArray(mechbcnodefacetbl, vectbl, {}, 'optpureproduct', true);
            
            linformvals = globlinformvals(mechbcinds, :);

            map = TensorMap();
            map.fromTbl = globmechbcnodefacevectbl;
            map.toTbl = mechbcnodefacevectbl;
            map.mergefds = {'bcinds', 'vec', 'nodes', 'faces'};
            map = map.setup();
            
            linform = map.eval(globlinform);
            linform = reshape(linform, dim, [])';
            
            % we need to remove 'bcinds' field because otherwise they will interfer with later assignment.
            mechbcnodefacetbl2 = replacefield(mechbcnodefacetbl, {{'bcinds', ''}});
            mechbc = struct('bcnodefacetbl', mechbcnodefacetbl2, ...
                            'linform'      , linform      , ...
                            'linformvals'  , linformvals);
        else

            mechbcnodefacetbl2 = replacefield(mechbcnodefacetbl, {{'bcinds', ''}});
            mechbc = struct('bcnodefacetbl', mechbcnodefacetbl2, ...
                            'linform'      , []      , ...
                            'linformvals'  , []);
        end

        % We get the part of external and volumetric force that are active in the block
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
        
        % We collect the degrees of freedom in the current block that belongs to the boundary for the fluid part.
        
        [fluidbcnodefacetbl, indstruct] = crossIndexArray(globfluidbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        globfluidbcnodeface_from_fluidbcnodeface = indstruct{1}.inds;
        
        bcterm_exists = true;
        if fluidbcnodefacetbl.num == 0
            bcterm_exists = false;
        end
        
        if bcterm_exists
            
            map = TensorMap();
            map.fromTbl  = globfluidbcnodefacetbl;
            map.toTbl    = fluidbcnodefacetbl;
            map.mergefds = {'faces', 'nodes'};
            map = map.setup();

            globfluidbcnodeface_from_fluidbcnodeface = map.getDispatchInd();
            
            fluidbcvals = map.eval(globfluidbcvals);

            clear bcdirichlet;
            bcdirichlet.bcnodefacetbl = fluidbcnodefacetbl;
            bcdirichlet.bcvals = fluidbcvals;
            
        else
            
            bcdirichlet = [];
            
        end
        
        % We get the part of external and volumetric sources that are active in the block
        map = TensorMap();
        map.fromTbl  = globnodefacetbl;
        map.toTbl    = nodefacetbl;
        map.mergefds = {'nodes', 'faces'};
        map = map.setup();
        
        extflux = map.eval(globextflux);
        
        % Assemble mechanical part
        
        map = TensorMap();
        map.fromTbl  = globfacetbl;
        map.toTbl    = facetbl;
        map.mergefds = {'faces'};
        map = map.setup();
        
        nnpf = map.eval(nnodesperface);
        
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero);
        output = coreMpsaAssembly(G, C, mechbc, nnpf, tbls, mappings, opts, 'useVirtual', useVirtual);

        mechmat    = output.matrices;
        mechbcvals = output.bcvals;
        
        % Assemble fluid part

        dooptimize = useVirtual;
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero, ...
                      'dooptimize', dooptimize);
        output = coreMpfaAssembly(G, K, bcdirichlet, tbls, mappings, opts, 'useVirtual', useVirtual);

        fluidmat    = output.matrices;
        fluidbcvals = output.bcvals;
        
        atbls = {nodefacevectbl, cellvectbl, nodefacetbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};
        locA = cell(nA, 1);
        for i = 1 : nA
            locA{i} = cell(4, 1);
            ni = atbls{i}.num;
            for j = 1 : nA
                nj = atbls{j}.num;
                locA{i}{j} = sparse(ni, nj);
            end
        end


        % Assemble coupling terms (finite volume and consistent divergence operators)
        
        map = TensorMap();
        map.fromTbl = globcelltbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        map = map.setup();
        
        alpha = map.eval(globalpha);
        rho   = map.eval(globrho);
        nnpc  = map.eval(nnodespercell);
        coupassembly = assembleCouplingTerms(G, eta, alpha, nnpc, tbls, mappings, 'useVirtual', useVirtual);
        
        % Recover terms from mpsa assembly
        
        invA11 = mechmat.invA11;
        locA{1}{1} = mechmat.A11;
        locA{1}{2} = mechmat.A12;
        locA{2}{1} = mechmat.A21;
        locA{2}{2} = mechmat.A22;
        locA{1}{5} = -mechmat.D;
        locA{5}{1} = -locA{1}{5}';
        
        % Recover terms from mpfa assembly
        
        invA33 = fluidmat.invA11;
        locA{3}{3} = fluidmat.A11;
        locA{3}{4} = fluidmat.A12;
        locA{4}{3} = fluidmat.A21;
        locA{4}{4} = fluidmat.A22;
        locA{3}{6} = -fluidmat.D;
        locA{6}{3} = -locA{3}{6}';

        % Recover the coupling terms
        
        locA{1}{4} = coupassembly.divfv;
        locA{1}{4} = -locA{1}{4}'; % We use the gradient which is the transpose of minus div
        locA{4}{1} = coupassembly.divconsnf;
        locA{4}{2} = coupassembly.divconsc;
    
        % We add the diagonal term for the mass conservation equation
        prod = TensorProd();
        prod.tbl1 = celltbl;
        prod.tbl2 = globcelltbl;
        prod.tbl3 = celltbl;
        prod.mergefds = {'cells'};
        prod = prod.setup();
        
        rho = prod.eval(rho, G.cells.volumes);
        
        % This matrix could be easily assembled directly (not using tensor assembly)
        celltbl = tbls.celltbl;
        prod = TensorProd();
        prod.tbl1 = celltbl;
        prod.tbl2 = celltbl;
        prod.tbl3 = celltbl;
        prod.mergefds = {'cells'};
        prod = prod.setup();
        
        A44b_T = SparseTensor();
        A44b_T = A44b_T.setFromTensorProd(rho, prod);
        locA{4}{4} = locA{4}{4} + A44b_T.getMatrix();
        
        % boundary conditions for the full system The volumetric force term (for mechanic) and source term (for fluid) is
        % added at the end, after the loop on the blocks.
        fullrhs = cell(6, 1);
        fullrhs{1} = extforce;
        fullrhs{3} = extflux;
        fullrhs{5} = mechbcvals;
        fullrhs{6} = fluidbcvals;

        % We proceed with the local reduction
        %
        %          | locB{1}{1}  locB{1}{2}  locB{1}{3}             |
        %  locB =  | locB{2}{1}  locB{2}{2}  locB{2}{3}  locB{2}{4} |
        %          | locB{3}{1}  locB{3}{2}  locB{3}{3}             |
        %          |             locB{4}{2}              locB{4}{4} |
        
        
        btbls  = {cellvectbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};        
        bnums  = cell(4, 1);
        locB   = cell(4, 4);
        locrhs = cell(4, 1);
        
        for i = 1 : 4
            locB{i} = cell(4, 1);
            ni = btbls{i}.num;
            bnums{i} = ni;
            locrhs{i} = zeros(ni, 1);
            for j = 1 : 4;
                ni = btbls{j}.num;
                locB{i}{j} = sparse(ni, nj);
            end
        end
    
        % 1. row : Momentum equation
        A21invA11 = locA{2}{1}*invA11;
        locB{1}{1} = -A21invA11*locA{1}{2} +  locA{2}{2};
        locB{1}{2} = -A21invA11*locA{1}{4};
        locB{1}{3} = -A21invA11*locA{1}{5};

        % 2. row : Fluid mass conservation
        A41invA11 = locA{4}{1}*invA11;
        A43invA33 = locA{4}{3}*invA33;
        locB{2}{1} = -A41invA11*locA{1}{2} + locA{4}{2};
        locB{2}{2} = -A41invA11*locA{1}{4} - A43invA33*locA{3}{4} + locA{4}{4};
        locB{2}{3} = -A41invA11*locA{1}{5};
        locB{2}{4} = -A43invA33*locA{3}{6};

        % 3. row : Mechanic BC
        A51invA11 = locA{5}{1}*invA11;
        locB{3}{1} = -A51invA11*locA{1}{2};
        locB{3}{2} = -A51invA11*locA{1}{4};
        locB{3}{3} = -A51invA11*locA{1}{5};

        % 4. row : Fluid BC
        A63invA33 = locA{6}{3}*invA33;
        locB{4}{2} = -A63invA33*locA{3}{4};
        locB{4}{4} = -A63invA33*locA{3}{6};

        % Assembly of right hand side
        f = fullrhs; % shortcut
        locrhs{1} = - A21invA11*f{1};
        locrhs{2} = - A41invA11*f{1} - A43invA33*f{3};
        locrhs{3} = f{5} - A51invA11*f{1};
        locrhs{4} = f{6} - A63invA33*f{3};
        
        % We store in the global matrices wthe values that have been computed at the block level
        
        btbls = {cellvectbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};        
        
        % Setup the index mappings (l2ginds) between local and global (could have been done more efficiently)
        fds = {{'cells', 'vec'}, ...
               {'cells'}, ...
               {'faces', 'nodes', 'bcinds'}, ...
               {'faces', 'nodes'}};
        l2ginds = cell(nB, 1);
        bnums = cell(nB, 1);
        
        for iInd = 1 : nB

            if ~isempty(btbls{iInd})
                
                map = TensorMap();
                map.fromTbl  = bglobtbls{iInd};
                map.toTbl    = btbls{iInd};
                map.mergefds = fds{iInd};

                if useVirtual
                    % Recall:
                    % bglobtbls = {globcellvectbl, globcelltbl, globmechbcnodefacetbl, globfluidbcnodefacetbl};    
                    % btbls     = {cellvectbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};
                    switch iInd
                      case 1
                        map.pivottbl = cellvectbl;
                        [vec, i] = ind2sub([vectbl.num, celltbl.num], (1 : cellvectbl.num)');
                        map.dispind1 = sub2ind([vectbl.num, globcelltbl.num], vec, globcell_from_cell(i));
                        map.dispind2 = (1 : cellvectbl.num)';
                      case 2
                        map.pivottbl = celltbl;
                        map.dispind1 = globcell_from_cell;
                        map.dispind2 = (1 : celltbl.num)';                    
                      case 3
                        map.pivottbl = mechbcnodefacetbl;
                        map.dispind1 = globmechbcnodeface_from_mechbcnodeface;
                        map.dispind2 = (1 : mechbcnodefacetbl.num)';
                      case 4
                        map.pivottbl = fluidbcnodefacetbl;
                        map.dispind1 = globfluidbcnodeface_from_fluidbcnodeface;
                        map.dispind2 = (1 : fluidbcnodefacetbl.num)';
                      otherwise
                        error('out of bound');
                    end
                    map.issetup = true;
                else
                    map = map.setup();
                end
                
                l2ginds{iInd} = map.getDispatchInd();
                
            else
                l2ginds{iInd} = [];
            end
        end
        
        for i = 1 : nB
            for j = 1 : nB
                [indi, indj, v] = find(locB{i}{j});
                indi = l2ginds{i}(indi);
                indj = l2ginds{j}(indj);
                ni = bglobnums{i};
                nj = bglobnums{j};
                insB = sparse(indi, indj, v, ni, nj); 
                B{i}{j} = B{i}{j} + insB;
            end
            [indi, ~, v] = find(locrhs{i});
            indi = l2ginds{i}(indi);
            ni = bglobnums{i};
            insrhs = sparse(indi, 1, v, ni, 1);
            rhs{i} = rhs{i} + insrhs;
        end
        
        if opt.assemblyMatrices
            % setup the index mappings (l2ginds)
            fds = {{'nodes', 'faces', 'vec'}   , ...
                   {'cells', 'vec'}            , ...
                   {'nodes', 'faces'}          , ...
                   {'cells'}                   , ...
                   {'faces', 'nodes', 'bcinds'}, ...
                   {'faces', 'nodes'}};
            
            al2ginds = cell(nA, 1);
            for i = 1 : nA
                map = TensorMap();
                map.fromTbl  = aglobtbls{i};
                map.toTbl    = atbls{i};
                map.mergefds = fds{i};
                
                al2ginds{i} = map.getDispatchInd();
            end
            
            for i = 1 : nA
                for j = 1 : nA
                    [indi, indj, v] = find(locA{i}{j});
                    indi = al2ginds{i}(indi);
                    indj = al2ginds{j}(indj);
                    ni = aglobtbls{i}.num;
                    nj = aglobtbls{j}.num;
                    insA{i}{j} = sparse(indi, indj, v, ni, nj); 
                    A{i}{j} = A{i}{j} + insA{i}{j};
                end
            end
        end
        
    end
    
    % We concatenate the matrices
    for i = 1 : nB
        B{i} = horzcat(B{i}{:});
    end
    B = vertcat(B{:});
    rhs{1} = rhs{1} + globforce;
    rhs{2} = rhs{2} + globsrc;
    
    rhs = vertcat(rhs{:});

    assembly = struct('B'  , B, ...
                      'rhs', rhs);
    
    if opt.assemblyMatrices
        assembly.A = A;
    end
    

end

