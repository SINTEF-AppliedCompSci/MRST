function output = runBiotConvSim(G, params, varargin)
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


    opt = struct('verbose'   , mrstVerbose, ...
                 'blocksize' , []         , ...
                 'useVirtual', false      , ...
                 'bcetazero' , false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;

    force_fun = params.force_fun;
    stress_fun = params.stress_fun;
    src_fun = params.src_fun;
    u_fun = params.u_fun;
    p_fun = params.p_fun;

    eta = params.eta; % continuity point

    lambda = params.lambda;
    mu     = params.mu;
    K      = params.K; % isotropic, diagonal coefficient
    alpha  = params.alpha;
    rho    = params.rho;

    [tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);

    coltbl         = tbls.coltbl;
    colrowtbl      = tbls.colrowtbl;
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;

    Nc = G.cells.num; 
    Nf = G.faces.num; 
    Nn = G.nodes.num; 
    Nd = G.griddim; 
    vdim = Nd*(Nd + 1)/2; % voigt dimension

    % prepare input for analytical functions
    for idim = 1 : Nd
        cc{idim} = G.cells.centroids(:, idim);
    end

    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);

    % Detect faces at each side
    xmin = min(G.faces.centroids(bcfaces, 1));
    isxmin = G.faces.centroids(bcfaces, 1) < (xmin + eps);
    xmax = max(G.faces.centroids(bcfaces, 1));
    isxmax = G.faces.centroids(bcfaces, 1) > (xmax - eps);
    ymin = min(G.faces.centroids(bcfaces, 2));
    isymin = G.faces.centroids(bcfaces, 2) < (ymin + eps);
    ymax = max(G.faces.centroids(bcfaces, 2));
    isymax = G.faces.centroids(bcfaces, 2) > (ymax - eps);
    if Nd == 3
        xmin = min(G.faces.centroids(bcfaces, 3));
        isxmin = G.faces.centroids(bcfaces, 3) < (xmin + eps);
        xmax = max(G.faces.centroids(bcfaces, 3));
        isxmax = G.faces.centroids(bcfaces, 3) > (xmax - eps);
    end

    isdir = (isxmin);
    dirfaces  = bcfaces(isdir);  % Dirichlet faces
    neumfaces = bcfaces(~isdir); % Neumann faces

    bcdirfacetbl.faces = dirfaces;
    bcdirfacetbl = IndexArray(bcdirfacetbl);
    bcdirnodefacetbl = crossIndexArray(bcdirfacetbl, nodefacetbl, {'faces'});
    bcdirnodefacecoltbl = crossIndexArray(bcdirnodefacetbl, coltbl, {}, 'optpureproduct', true);

    [~, nodefacecents] = computeNodeFaceCentroids(G, eta, tbls, 'bcetazero', opt.bcetazero);

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcdirnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();

    bcdirnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacecoltbl:
    bcdirnodefacecents = reshape(bcdirnodefacecents, Nd, [])';
    bcnum = bcdirnodefacetbl.num;

    % Prepare input for analytical functions
    for idim = 1 : Nd
        bnfc{idim} = bcdirnodefacecents(:, idim);
    end

    % Compute boundary conditions for the mechanical part (Dirichlet and Neumann parts)

    % 1) Compute Dirichlet part of mechanical bc
    for idim = 1 : Nd
        linform = zeros(bcnum, Nd);
        linform(:, idim) = 1;
        linforms{idim} = linform;
        linformvals{idim} = u_fun{idim}(bnfc{:});
    end

    bcdirfaces = bcdirnodefacetbl.get('faces');
    bcdirnodes = bcdirnodefacetbl.get('nodes');
    extbcdirnodefacetbl.faces = repmat(bcdirfaces, Nd, 1);
    extbcdirnodefacetbl.nodes = repmat(bcdirnodes, Nd, 1);
    extbcdirnodefacetbl = IndexArray(extbcdirnodefacetbl);

    bc.bcnodefacetbl = extbcdirnodefacetbl;
    bc.linform = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});
    clear extbcnodefacetbl linforms linformvals

    % 2) Compute Neumann part of mechanical bc

    bcneumfacetbl.faces = neumfaces;
    bcneumfacetbl = IndexArray(bcneumfacetbl);
    bcneumnodefacetbl = crossIndexArray(bcneumfacetbl, nodefacetbl, {'faces'});
    bcneumnodefacecoltbl = crossIndexArray(bcneumnodefacetbl, coltbl, {}, 'optpureproduct', true);
    voigttbl.voigtdim = (1 : vdim)';
    voigttbl = IndexArray(voigttbl);
    bcneumnodefacevoigttbl = crossIndexArray(bcneumnodefacetbl, voigttbl, {}, 'optpureproduct', true);    
    colrowvoigttbl = colrowtbl;
    switch Nd
      case 2
        voigtind = [1; 3; 3; 2];
      case 3
        voigtind = [1; 6; 5; 6; 2; 4; 5; 4; 3];
      otherwise
        error('d not recognized');
    end
    colrowvoigttbl = colrowvoigttbl.addInd('voigtdim', voigtind);

    bcneumnodefacevoigttbl = crossIndexArray(bcneumnodefacetbl, voigttbl, {}, 'optpureproduct', true);
    bcneumnodefacecolrowvoigttbl = crossIndexArray(bcneumnodefacetbl, colrowvoigttbl, {}, 'optpureproduct', true);

    % Evaluate analytical stress at bcneumnodefacecents (continuity points at boundary)

    % We get the continuity points where the stress should be evaluated, they belong to bcnemnodefacecoltbl
    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcneumnodefacecoltbl;
    map.mergefds = {'nodes',  'faces', 'coldim'};
    map = map.setup();

    bcneumnodefacecents = map.eval(nodefacecents);
    bcneumnodefacecents = reshape(bcneumnodefacecents, Nd, [])';
    for idim = 1 : Nd
        bnfc{idim} = bcneumnodefacecents(:, idim);
    end
    for idim = 1 : voigttbl.num
        stress{idim} = stress_fun{idim}(bnfc{:});
    end
    pneum = p_fun(bnfc{:});

    % Reformat stress in bcneumnodefacevoigttbl
    stress = horzcat(stress{:});
    stress = reshape(stress', [], 1);

    % Map stress to bcneumnodefacecolrowvoigttbl
    map = TensorMap();
    map.fromTbl = bcneumnodefacevoigttbl;
    map.toTbl = bcneumnodefacecolrowvoigttbl;
    map.mergefds = {'nodes', 'faces', 'voigtdim'};
    map = map.setup();

    stress = map.eval(stress);

    % Fetch the normals in bcneumnodefacecoltbl
    cellnodefacetbl = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    % facetNormals is in cellnodefacecoltbl;
    facetNormals =  computeFacetNormals(G, cellnodefacetbl);

    map = TensorMap();
    map.fromTbl = cellnodefacecoltbl;
    map.toTbl = bcneumnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    % now facetNormals is in bcneumnodefacecoltbl
    facetNormals = map.eval(facetNormals);

    % We multiplty stress with facetNormals
    prod = TensorProd();
    prod.tbl1 = bcneumnodefacecolrowvoigttbl;
    prod.tbl2 = bcneumnodefacecoltbl;
    prod.tbl3 = bcneumnodefacecoltbl;
    prod.replacefds1 = {{'coldim', 'rowdim', 'interchange'}};
    prod.replacefds2 = {{'coldim', 'rowdim'}};
    prod.mergefds = {'nodes', 'faces'};
    prod.reducefds = {'rowdim'};
    prod = prod.setup();

    extforce = prod.eval(stress, facetNormals);

    % We add the part due to the Biot pressure
    apneum = - alpha*pneum;

    prod = TensorProd();
    prod.tbl1 = bcneumnodefacetbl;
    prod.tbl2 = bcneumnodefacecoltbl;
    prod.tbl3 = bcneumnodefacecoltbl;
    prod.mergefds = {'nodes', 'faces'};
    prod = prod.setup();

    extforce = extforce + prod.eval(apneum, facetNormals);

    % the format of extforce expected in assembly is in nodefacecoltbl
    map = TensorMap();
    map.fromTbl = bcneumnodefacecoltbl;
    map.toTbl = nodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();

    extforce = map.eval(extforce);

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
    loadstruct.extforce = extforce;
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
    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', true);

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();

    bcnodefacecents = map.eval(nodefacecents);

    % Here, we assume a given structure of bcnodefacecoltbl:
    bcnodefacecents = reshape(bcnodefacecents, Nd, [])';
    % Prepare input for analytical functions
    for idim = 1 : Nd
        bnfc{idim} = bcnodefacecents(:, idim);
    end
    bcpvals = p_fun(bnfc{:});
    bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                         'bcvals'       , bcpvals);
    bcstruct = struct('bcdirichlet', bcdirichlet, ...
                      'bcneumann'  , bcneumann);
    src = src_fun(cc{:});
    src = G.cells.volumes.*src;
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

    K = K*eye(Nd);
    K = reshape(K, [], 1);
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
        options = {'verbose'  , true         , ...
                   'bcetazero', opt.bcetazero, ...
                   'blocksize', opt.blocksize};
        assembly = blockAssembleBiot(G, props, drivingforces, eta, tbls, mappings, options{:});
    else
        assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, 'bcetazero', opt.bcetazero);
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
