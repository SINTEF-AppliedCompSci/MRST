function assembly = blockAssembleBiot(G, props, drivingforces, eta, globtbls, globmappings, varargin)
    
    opt = struct('verbose'         , mrstVerbose, ...
                 'assemblyMatrices', false      , ...
                 'addAdOperators'  , false      , ...
                 'blocksize'       , []         , ...
                 'bcetazero'       , true       , ...
                 'useVirtual'      , true       , ...
                 'extraoutput'     , false);
        
    opt = merge_options(opt, varargin{:});
    
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
    %       | displacement_nfc (node face dofs, belongs to nodefacecoltbl)         |
    %       | displacement_c   (cell dofs, belongs to cellcoltbl)                  |
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
    blocksize = opt.blocksize;
    
    nn = G.nodes.num;
    nblocks = floor(nn/blocksize);
    blocksizes = repmat(blocksize, nblocks, 1);
    if nn > nblocks*blocksize
        blocksizes = [blocksizes; nn - nblocks*blocksize];
    end
    nblocks = numel(blocksizes);
    blockinds = cumsum([1; blocksizes]);

    coltbl         = globtbls.coltbl;
    colrowtbl      = globtbls.colrowtbl;
    col2row2tbl    = globtbls.col2row2tbl;
    
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
    
    dim = coltbl.num;
    
    
    mechprops  = props.mechprops;
    fluidprops = props.fluidprops;
    coupprops  = props.coupprops;
   
    globC     = setupStiffnessTensor(mechprops, globtbls);
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
    globmechbcnodefacecoltbl = crossIndexArray(globmechbcnodefacetbl, coltbl, {}, ...
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
    %       | displacement at cell          |    in globcellcoltbl
    %  u =  | pressure at cell              |    in globcelltbl
    %       | lagrange multiplier mechanics |    in globmechbcnodefacetbl
    %       | lagrange multiplier fluid     |    in globfluidbcnodefacetbl
    

    bglobtbls = {globcellcoltbl, globcelltbl, globmechbcnodefacetbl, globfluidbcnodefacetbl};
    
    B   = cell(4, 1);
    rhs = cell(4, 1);
    bglobnums = cell(4, 1);
    
    for i = 1 : 4
        B{i} = cell(4, 1);
        ni = bglobtbls{i}.num;
        bglobnums{i} = ni;
        rhs{i} = zeros(ni, 1);
        for j = 1 : 4;
            nj = bglobtbls{j}.num;
            B{i}{j} = sparse(ni, nj);
        end
    end
    
    for iblock = 1 : nblocks

        %% Construction of tensor g (as defined in paper eq 4.1.2)
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
    
        % Obtain stiffness values for the block
        map = TensorMap();
        map.fromTbl = globcellcol2row2tbl;
        map.toTbl = cellcol2row2tbl;
        map.mergefds = {'cells', 'coldim1', 'coldim2', 'rowdim1', 'rowdim2'};
        
        map.pivottbl = cellcol2row2tbl;
        c_num     = celltbl.num;             % shortcut
        gc_num    = globcellcol2row2tbl.num; % shortcut
        cc2r2_num = cellcol2row2tbl.num;     % shortcut
        c2r2_num  = col2row2tbl.num;         % shortcut
        [c2r2, i] = ind2sub([c2r2_num, c_num], (1 : cc2r2_num)');
        map.dispind1 = sub2ind([c2r2_num, gc_num], c2r2, globcell_from_cell(i));
        map.dispind2 = (1 : cc2r2_num)';
        map.issetup = true;
        
        C = map.eval(globC);

        % Obtain permeability values for the block
        map = TensorMap();
        map.fromTbl = globcellcolrowtbl;
        map.toTbl = cellcolrowtbl;
        map.mergefds = {'cells', 'coldim', 'rowdim'};
        
        map.pivottbl = cellcolrowtbl;
        ccr_num = cellcolrowtbl.num; % shortcut
        cr_num = colrowtbl.num;      % shortcut
        [cr, i] = ind2sub([cr_num, c_num], (1 : ccr_num)');
        map.dispind1 = sub2ind([cr_num, gc_num], cr, globcell_from_cell(i));
        map.dispind2 = (1 : ccr_num)';
        map.issetup = true;
        
        K = map.eval(globK);
        

        %% Assembly mechanical part
        
        % We collect the degrees of freedom in the current block that belongs to the boundary.
        mechbcnodefacetbl = crossIndexArray(globmechbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        
        mechbcterm_exists = true;
        if mechbcnodefacetbl.num == 0
            mechbcterm_exists = false;
        end
                
        if mechbcterm_exists
            
            mechbcinds = mechbcnodefacetbl.get('bcinds');
            mechbcnodefacecoltbl = crossIndexArray(mechbcnodefacetbl, coltbl, {}, 'optpureproduct', true);
            
            linformvals = globlinformvals(mechbcinds, :);

            map = TensorMap();
            map.fromTbl = globmechbcnodefacecoltbl;
            map.toTbl = mechbcnodefacecoltbl;
            map.mergefds = {'bcinds', 'coldim', 'nodes', 'faces'};
            map = map.setup();
            
            linform = map.eval(globlinform);
            linform = reshape(linform, dim, [])';
            
            % we need to remove 'bcinds' field because otherwise they will interfer with later assignment.
            mechbcnodefacetbl2 = replacefield(mechbcnodefacetbl, {{'bcinds', ''}});
            mechbc = struct('bcnodefacetbl', mechbcnodefacetbl2, ...
                            'linform'      , linform      , ...
                            'linformvals'  , linformvals);
        end
        
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero);
        [mechmat, mechbcvals] = coreMpsaAssembly(G, C, mechbc, tbls, mappings, opts);
        
        invA11 = mechmat.invA11;
        A11 = mechmat.A11;
        A12 = mechmat.A12;
        A21 = mechmat.A21;
        A22 = mechmat.A22;
        A15 = -mechmat.D;
        A51 = -A15';

        % We get the part of external and volumetric force that are active in the block
        map = TensorMap();
        map.fromTbl = globnodefacecoltbl;
        map.toTbl = nodefacecoltbl;
        map.mergefds = {'nodes', 'faces', 'coldim'};
        map = map.setup();

        extforce = map.eval(globextforce);
        
        map = TensorMap();
        map.fromTbl = globcellcoltbl;
        map.toTbl = cellcoltbl;
        map.mergefds = {'cells', 'coldim'};
        map = map.setup();

        force = map.eval(globforce);
        
        %% Assembly of fluid part
        
        % We collect the degrees of freedom in the current block that belongs to the boundary.
        
        fluidbcnodefacetbl = crossIndexArray(globfluidbcnodefacetbl, nodefacetbl, {'nodes', 'faces'});
        
        bcterm_exists = true;
        if fluidbcnodefacetbl.num == 0
            bcterm_exists = false;
        end
        
        if bcterm_exists
            map = TensorMap();
            map.fromTbl = globfluidbcnodefacetbl;
            map.toTbl = fluidbcnodefacetbl;
            map.mergefds = {'faces', 'nodes'};
            map = map.setup();
            
            fluidbcvals = map.eval(globfluidbcvals);

            clear bcdirichlet;
            bcdirichlet.bcnodefacetbl = fluidbcnodefacetbl;
            bcdirichlet.bcvals = fluidbcvals;
        else
            bcdirichlet = [];
        end
        
        dooptimize = useVirtual;
        opts = struct('eta', eta, ...
                      'bcetazero', opt.bcetazero, ...
                      'dooptimize', dooptimize);
        [fluidmat, fluidbcvals, extra] = coreMpfaAssembly(G, K, bcdirichlet, tbls, mappings, opts);
        
        invA33 = fluidmat.invA11;
        A33 = fluidmat.A11;
        A34 = fluidmat.A12;
        A43 = fluidmat.A21;
        A44 = fluidmat.A22;
        A36 = -fluidmat.D;
        A63 = -A36';
        
        % We get the part of external and volumetric sources that are active in the block
        map = TensorMap();
        map.fromTbl = globnodefacetbl;
        map.toTbl = nodefacetbl;
        map.mergefds = {'nodes', 'faces'};
        map = map.setup();
        
        extflux = map.eval(globextflux);
        
        map = TensorMap();
        map.fromTbl = globcelltbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        map = map.setup();

        src = map.eval(globsrc);

        %% Assemble coupling terms (finite volume and consistent divergence operators)
        
        map = TensorMap();
        map.fromTbl = globcelltbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        map = map.setup();
        
        alpha = map.eval(globalpha);
        rho = map.eval(globrho);
        
        coupassembly = assembleCouplingTerms(G, eta, alpha, tbls, mappings);
        
        % Recover the coupling terms
        A14 = coupassembly.divfv;
        A14 = -A14'; % We use the gradient which is the transpose of minus div
        A41 = coupassembly.divconsnf;
        A42 = coupassembly.divconsc;
    
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
        A44 = A44 + A44b_T.getMatrix();
        
        % boundary conditions for the full system
        fullrhs = cell(6, 1);
        fullrhs{1} = extforce;
        fullrhs{2} = force;
        fullrhs{3} = extflux;
        fullrhs{4} = src;
        fullrhs{5} = mechbcvals;
        fullrhs{6} = fluidbcvals;

        %% We proceed with the local reduction
        
        %          | locB{1}{1}  locB{1}{2}  locB{1}{3}             |
        %  locB =  | locB{2}{1}  locB{2}{2}  locB{2}{3}  locB{2}{4} |
        %          | locB{3}{1}  locB{3}{2}  locB{3}{3}             |
        %          |             locB{4}{2}              locB{4}{4} |
        
        
        btbls  = {cellcoltbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};        
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
        A21invA11 = A21*invA11;
        locB{1}{1} = -A21invA11*A12 +  A22;
        locB{1}{2} = -A21invA11*A14;
        locB{1}{3} = -A21invA11*A15;

        % 2. row : Fluid mass conservation
        A41invA11 = A41*invA11;
        A43invA33 = A43*invA33;
        locB{2}{1} = -A41invA11*A12 + A42;
        locB{2}{2} = -A41invA11*A14 - A43invA33*A34 + A44;
        locB{2}{3} = -A41invA11*A15;
        locB{2}{4} = -A43invA33*A36;

        % 3. row : Mechanic BC
        A51invA11 = A51*invA11;
        locB{3}{1} = -A51invA11*A12;
        locB{3}{2} = -A51invA11*A14;
        locB{3}{3} = -A51invA11*A15;

        % 4. row : Fluid BC
        A63invA33 = A63*invA33;
        locB{4}{2} = -A63invA33*A34;
        locB{4}{4} = -A63invA33*A36;

        % Assembly of right hand side
        f = fullrhs; % shortcut
        locrhs{1} = f{2} - A21invA11*f{1};
        locrhs{2} = f{4} - A41invA11*f{1} - A43invA33*f{3};
        locrhs{3} = f{5} - A51invA11*f{1};
        locrhs{4} = f{6} - A63invA33*f{3};
        
        % We store in the global matrices wthe values that have been computed at the block level
        
        btbls = {cellcoltbl, celltbl, mechbcnodefacetbl, fluidbcnodefacetbl};        
        
        % setup the index mappings (l2ginds)
        fds = {{'cells', 'coldim'}, {'cells'}, {'faces', 'nodes', 'bcinds'}, {'faces', 'nodes'}};
        l2ginds = cell(4, 1);
        bnums = cell(4, 1);
        for i = 1 : 4
            bnums{i} = btbls{i}.num;
            map = TensorMap();
            map.fromTbl = bglobtbls{i};
            map.toTbl = btbls{i};
            map.mergefds = fds{i};
            
            l2ginds{i} = map.getDispatchInd();
        end
        
        for i = 1 : 4
            for j = 1 : 4
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
    end
    
    % We concatenate the matrices
    for i = 1 : 4
        B{i} = horzcat(B{i}{:});
    end
    B = vertcat(B{:});
    rhs = vertcat(rhs{:});

    assembly = struct('B'  , B, ...
                      'rhs', rhs);
    
    if opt.assemblyMatrices
        fullsystem.A = A;
        fullsystem.rhs = cat(fullrhs{:});
        assembly.fullsystem = fullsystem;
    end    
    
    if opt.addAdOperators

        fluxop = fluidassembly.adoperators.fluxop;
        
        % Setup face node dislpacement operator
        fndisp{1} = -invA11*A12;
        fndisp{2} = -invA11*A14;
        fndisp{3} = -invA11*A15;
        
        facenodedispop = @(u, p, lm, extforce) facenodedispopFunc(u, p, lm, extforce, fndisp);
        
        % Setup stress operator
        aver = cellAverageOperator(tbls, mappings);
        stress{1} = C1;
        stress{2} = C2;
        stressop = @(unf, uc) stressopFunc(unf, uc, stress, aver);

        % Setup divKgrad operator
        divKgrad{1} = - A43invA33*A34 + A44;
        divKgrad{2} = - A43invA33*A36;
        divKgradrhs = f{4} - A43invA33*f{3};
        
        divKgradop = @(p, lf) divKgradopFunc(p, lf, divKgrad, divKgradrhs);

        % Setup consistent divergence operator (for displacement, includes value of Biot coefficient alpha)
        % The divergence is volume weighted
        divu{1} = - A41invA11*A12 + A42;
        divu{2} = - A41invA11*A14;
        divu{3} = - A41invA11*A15;
        divu{4} = A41invA11;
        
        divuop = @(u, p, lm, extforce) divuopFunc(u, p, lm, extforce, divu);

        % Setup momentum balance operator 
        moment{1} = B11;
        moment{2} = B12;
        moment{3} = B13;
        % We have : right-hand side for momentum equation = f{2} - A21invA11*f{1}. Hence, we set
        moment{4} = A21invA11;
        momentrhs = f{2};
        
        momentop = @(u, p, lm, extforce) momentopFunc(u, p, lm, extforce, moment, momentrhs);
        
        % Setup dirichlet boundary operator for mechanics
        mechdir{1} = B31;
        mechdir{2} = B32;
        mechdir{3} = B33;
        mechdirrhs = redrhs{3};
        
        mechDirichletop = @(u, p, lm) mechdiropFunc(u, p, lm, mechdir, mechdirrhs);
        
        % Setup dirichlet boundary operator for flow
        fluiddir{1} = B42;
        fluiddir{2} = B44;
        fluiddirrhs = redrhs{4};
        
        fluidDirichletop = @(p, lf) fluiddiropFunc(p, lf, fluiddir, fluiddirrhs);
        
        adoperators = struct('fluxop'          , fluxop          , ...
                             'facenodedispop'  , facenodedispop  , ...
                             'stressop'        , stressop        , ...
                             'divKgradop'      , divKgradop      , ...
                             'divuop'          , divuop          , ...
                             'momentop'        , momentop        , ...
                             'fluidDirichletop', fluidDirichletop, ...
                             'mechDirichletop' , mechDirichletop);

        assembly.adoperators = adoperators;
        
    end    
end


function fndisp = facenodedispopFunc(u, p, lm, extforce, fndisp)
    fndisp = fndisp{1}*u + fndisp{2}*p + fndisp{3}*lm + extforce;
end

function stress = stressopFunc(unf, uc, stress, aver)
    
    % get stress at each cell-node region (corner)
    stress = stress{1}*unf + stress{2}*uc;
    stress = aver*stress;
    
end

function divKgrad = divKgradopFunc(p, lf, divKgrad, divKgradrhs)
    divKgrad = divKgrad{1}*p + divKgrad{2}*lf - divKgradrhs;
end

function divu = divuopFunc(u, p, lm, extforce, divu)
    divu = divu{1}*u + divu{2}*p + divu{3}*lm + divu{4}*extforce;
end

function moment = momentopFunc(u, p, lm, extforce, moment, momentrhs)
    moment = moment{1}*u + moment{2}*p + moment{3}*lm + moment{4}*extforce - momentrhs;
end

function mechdir = mechdiropFunc(u, p, lm, mechdir, mechdirrhs)
    mechdir = mechdir{1}*u + mechdir{2}*p + mechdir{3}*lm - mechdirrhs;
end

function fluiddir = fluiddiropFunc(p, lf, fluiddir, fluiddirrhs)
    fluiddir = fluiddir{1}*p + fluiddir{2}*lf - fluiddirrhs;
end
% 