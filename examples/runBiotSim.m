function output = runBiotSim(G, params, varargin)
    
    opt = struct('verbose'   , mrstVerbose, ...
                 'blocksize' , []         , ...
                 'useVirtual', false      , ...
                 'bcetazero' , false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    force_fun = params.force_fun;
    src_fun = params.src_fun;
    u_fun = params.u_fun;
    p_fun = params.p_fun;
    
    eta = params.eta; % continuity point
    
    lambda = params.lambda;
    mu     = params.mu;
    K      = params.K; % isotropic, diagonal coefficient
    alpha  = params.alpha;
    rho    = params.rho;
   
    doVem = false;
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);
    coltbl = tbls.coltbl;
    nodefacetbl = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    
    Nc = G.cells.num; 
    Nf = G.faces.num; 
    Nn = G.nodes.num; 
    Nd = G.griddim; 
    
    % prepare input for analytical functions
    for idim = 1 : Nd
        cc{idim} = G.cells.centroids(:, idim);
    end
    
    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);
    
    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', ...
                                       true);
    clear bcfacetbl
    
    [~, nodefacecents] = computeNodeFaceCentroids(G, eta, tbls, 'bcetazero', opt.bcetazero);
    
    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    
    bcnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacecoltbl:
    bcnodefacecents = reshape(bcnodefacecents, Nd, [])';
    bcnum = bcnodefacetbl.num;
   
    dotest = false;
    if dotest
        % plot continuity points 
        figure
        hold on
        plotGrid(G)
        nodefacecents = reshape(nodefacecents, Nd, [])';
        if Nd == 2
            plot(nodefacecents(:, 1), nodefacecents(:, 2), '*');
        else
            plot3(nodefacecents(:, 1), nodefacecents(:, 2), nodefacecents(:, 3), '*');
        end
    end
    
    
    % Prepare input for analytical functions
    for idim = 1 : Nd
        bnfc{idim} = bcnodefacecents(:, idim);
    end
    
    % Compute boundary conditions for mechanical part
    for idim = 1 : Nd
        linform = zeros(bcnum, Nd);
        linform(:, idim) = 1;
        linforms{idim} = linform;
        linformvals{idim} = u_fun{idim}(bnfc{:});
    end
    
    bcfaces = bcnodefacetbl.get('faces');
    bcnodes = bcnodefacetbl.get('nodes');
    extbcnodefacetbl.faces = repmat(bcfaces, Nd, 1);
    extbcnodefacetbl.nodes = repmat(bcnodes, Nd, 1);
    extbcnodefacetbl = IndexArray(extbcnodefacetbl);
    
    bc.bcnodefacetbl = extbcnodefacetbl;
    bc.linform = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});
    clear extbcnodefacetbl linforms linformvals

    % Compute body force
    force = NaN(Nc, Nd);
    for idim = 1 : Nd
        force(:, idim) = force_fun{idim}(cc{:});
    end
    force = bsxfun(@times, G.cells.volumes, force);
    % Here, we assume we know the structure of cellcoltbl;
    force = reshape(force', [], 1);
    
    loadstruct.bc = bc;
    loadstruct.force = force;
    loadstruct.extforce = zeros(tbls.nodefacecoltbl.num, 1);
    clear bc force
    
    % Setup mechanical parameters
    mu = mu*ones(Nc, 1);
    lambda = lambda*ones(Nc, 1);
    mechprops = struct('mu'  , mu, ...
                       'lambda', lambda);
    
    % Setup boundary conditions for fluid part 
    % Neumann fluid bc
    bcneumann = [];
    % Dirichlet fluid bc
    bcpvals = p_fun(bnfc{:});
    bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                         'bcvals'       , bcpvals);
    bcstruct = struct('bcdirichlet', bcdirichlet, ...
                      'bcneumann'  , bcneumann);
    src = src_fun(cc{:});
    fluidforces = struct('bcstruct', bcstruct, ...
                         'src'     , src);
    
    % Setup fluid parameters
    cellcolrowtbl = tbls.cellcolrowtbl;
    colrowtbl = tbls.colrowtbl;

    map = TensorMap();
    map.fromTbl = colrowtbl;
    map.toTbl = cellcolrowtbl;
    map.mergefds = {'coldim', 'rowdim'};
    map = map.setup();

    K = [K; 0; 0; K];
    K = map.eval(K);
    fluidprops.K = K;

    % Setup coupling parameters
    alpha = alpha*ones(Nc, 1);
    rho = rho*ones(Nc, 1);
    coupprops = struct('alpha', alpha, ...
                       'rho', rho);
    
    props = struct('mechprops' , mechprops , ...
                   'fluidprops', fluidprops, ...
                   'coupprops' , coupprops);
    
    % setup driving forces
    drivingforces = struct('mechanics', loadstruct, ...
                           'fluid'    , fluidforces);
                           
  
    if ~isempty(opt.blocksize)
        error('not yet implemented');
    else
        assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings);
    end
    
    clear prop loadstruct
    
    B   = assembly.B;
    rhs = assembly.rhs;
    sol = B\rhs;

    % Displacement values at cell centers.
    cellcoltbl = tbls.cellcoltbl;
    ncc = cellcoltbl.num;

    u = sol(1 : ncc);
    u = reshape(u, Nd, [])';    
    
    % Pressure values at cell centers.
    celltbl = tbls.celltbl;
    nc = celltbl.num;

    p = sol(ncc + 1 : ncc + nc);

    output = struct('B'       , B       , ...
                    'assembly', assembly, ...
                    'rhs'     , rhs     , ...
                    'tbls'    , tbls    , ...
                    'u'       , u       , ...
                    'p'       , p);
    
end
