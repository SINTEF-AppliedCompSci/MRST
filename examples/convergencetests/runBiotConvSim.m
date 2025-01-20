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

    [tbls, mappings] = setupMpxaStandardTables(G, 'useVirtual', useVirtual);

    vectbl         = tbls.vectbl;
    vec12tbl       = tbls.vec12tbl;
    nodefacetbl    = tbls.nodefacetbl;
    nodefacevectbl = tbls.nodefacevectbl;

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

    bcdirfacetbl.faces  = dirfaces;
    bcdirfacetbl = IndexArray(bcdirfacetbl);
    [bcdirnodefacetbl, indstruct] = crossIndexArray(bcdirfacetbl, nodefacetbl, {'faces'});
    if useVirtual
        d_num   = vectbl.num;
        bcdirface_from_bcdirnodeface = indstruct{1}.inds;
        nodeface_from_bcdirnodeface  = indstruct{2}.inds;
    end
    bcdirnodefacevectbl = crossIndexArray(bcdirnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

    [~, nodefacecents] = computeNodeFaceCentroids2(G, eta, tbls              , ...
                                                   'bcetazero', opt.bcetazero, ...
                                                   'useVirtual', useVirtual);

    map = TensorMap();
    map.fromTbl  = nodefacevectbl;
    map.toTbl    = bcdirnodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual
        map.pivottbl = bcdirnodefacevectbl;

        N = bcdirnodefacetbl.num; 
        [vec, i] = ind2sub([d_num, N], (1 : bcdirnodefacevectbl.num)');

        map.dispind1 = sub2ind([d_num, nodefacetbl.num], vec , nodeface_from_bcdirnodeface(i));
        map.dispind2 = (1 : bcdirnodefacevectbl.num)';
        
        map.issetup = true;
        
    else
        map = map.setup();
    end

    
    bcdirnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacevectbl:
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
    [bcneumnodefacetbl, indstruct] = crossIndexArray(bcneumfacetbl, nodefacetbl, {'faces'});
    if useVirtual
        bcneumface_from_bcneumnodeface = indstruct{1}.inds;
        nodeface_from_bcneumnodeface   = indstruct{2}.inds;
    end
    bcneumnodefacevectbl = crossIndexArray(bcneumnodefacetbl, vectbl, {}, ...
                                           'optpureproduct' , true  , ...
                                           'virtual'        , useVirtual);
    voigttbl.voigtdim = (1 : vdim)';
    voigttbl = IndexArray(voigttbl);
    bcneumnodefacevoigttbl = crossIndexArray(bcneumnodefacetbl, voigttbl, {}, ...
                                             'optpureproduct', true         , ...
                                             'virtual', useVirtual);    
    vec12voigttbl = vec12tbl;
    switch Nd
      case 2
        voigt_from_vec12 = [1; 3; 3; 2];
      case 3
        voigt_from_vec12 = [1; 6; 5; 6; 2; 4; 5; 4; 3];
      otherwise
        error('d not recognized');
    end
    vec12voigttbl = vec12voigttbl.addInd('voigtdim', voigt_from_vec12);

    
    bcneumnodefacevec12voigttbl = crossIndexArray(bcneumnodefacetbl, vec12voigttbl, {}, ...
                                                  'optpureproduct', true              , ...
                                                  'virtual'    , useVirtual);

    % Evaluate analytical stress at bcneumnodefacecents (continuity points at boundary)

    % We get the continuity points where the stress should be evaluated, they belong to bcnemnodefacevectbl
    map = TensorMap();
    map.fromTbl  = nodefacevectbl;
    map.toTbl    = bcneumnodefacevectbl;
    map.mergefds = {'nodes',  'faces', 'vec'};

    if useVirtual
        
        map.pivottbl = bcneumnodefacevectbl;

        N = bcneumnodefacetbl.num; 
        [vec, i] = ind2sub([d_num, N], (1 : bcneumnodefacevectbl.num)');

        map.dispind1 = sub2ind([d_num, nodefacetbl.num], vec, nodeface_from_bcneumnodeface(i));
        map.dispind2 = (1 : bcneumnodefacevectbl.num)';
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end
    
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

    % Map stress to bcneumnodefacevec12voigttbl
    map = TensorMap();
    map.fromTbl  = bcneumnodefacevoigttbl;
    map.toTbl    = bcneumnodefacevec12voigttbl;
    map.mergefds = {'nodes', 'faces', 'voigtdim'};

    if useVirtual
        
        map.pivottbl = bcneumnodefacevec12voigttbl;

        N = bcneumnodefacetbl.num; 
        [v, i] = ind2sub([vec12tbl.num, N], (1 : bcneumnodefacevec12voigttbl.num)');

        map.dispind1 = sub2ind([voigttbl.num, bcneumnodefacetbl.num], voigt_from_vec12(v), i);
        map.dispind2 = (1 : bcneumnodefacevec12voigttbl.num)';
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end

    stress = map.eval(stress);

    % Fetch the normals in bcneumnodefacevectbl
    cellnodefacetbl = tbls.cellnodefacetbl;
    cellnodefacevectbl = tbls.cellnodefacevectbl;
    % facetNormals is in cellnodefacevectbl;
    facetNormals =  computeFacetNormals(G, cellnodefacetbl);

    map = TensorMap();
    map.fromTbl  = cellnodefacevectbl;
    map.toTbl    = bcneumnodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual

        [bcneumcellnodefacetbl, indstruct] = crossIndexArray(cellnodefacetbl, bcneumnodefacetbl, {'nodes', 'faces'});
        cellnodeface_from_bcneumcellnodefacetbl   = indstruct{1}.inds;
        bcneumnodeface_from_bcneumcellnodefacetbl = indstruct{2}.inds;
        bcneumcellnodefacevectbl = crossIndexArray(bcneumcellnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
        
        map.pivottbl = bcneumcellnodefacevectbl;

        N = bcneumcellnodefacetbl.num; 
        [vec, i] = ind2sub([d_num, N], (1 : bcneumcellnodefacevectbl.num)');

        map.dispind1 = sub2ind([d_num, cellnodefacetbl.num], vec, cellnodeface_from_bcneumcellnodefacetbl(i));
        map.dispind2 = sub2ind([d_num, bcneumnodefacetbl.num], vec, bcneumnodeface_from_bcneumcellnodefacetbl(i));
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end
    % now facetNormals is in bcneumnodefacevectbl
    facetNormals = map.eval(facetNormals);

    % We multiplty stress with facetNormals
    prod = TensorProd();
    prod.tbl1        = bcneumnodefacevec12voigttbl;
    prod.tbl2        = bcneumnodefacevectbl;
    prod.tbl3        = bcneumnodefacevectbl;
    prod.replacefds1 = {{'vec1', 'redvec'}, {'vec2', 'vec'}};
    prod.replacefds2 = {{'vec', 'redvec'}};
    prod.mergefds    = {'nodes', 'faces'};
    prod.reducefds   = {'redvec'};

    if useVirtual
        
        prod.pivottbl = bcneumnodefacevec12voigttbl;
        
        N = bcneumnodefacetbl.num;
        [vec2, vec1, i] = ind2sub([d_num, d_num, N], (1 : bcneumnodefacevec12voigttbl.num)');

        prod.dispind1 = (1 : bcneumnodefacevec12voigttbl.num)';
        N = bcneumnodefacetbl.num;
        prod.dispind2 = sub2ind([d_num, N], vec1, i);
        prod.dispind3 = sub2ind([d_num, N], vec2, i);
        
        prod.issetup = true;
        
    else
        
        prod = prod.setup();
        
    end
    
    extforce = prod.eval(stress, facetNormals);

    % We add the part due to the Biot pressure
    apneum = - alpha*pneum;

    prod = TensorProd();
    prod.tbl1 = bcneumnodefacetbl;
    prod.tbl2 = bcneumnodefacevectbl;
    prod.tbl3 = bcneumnodefacevectbl;
    prod.mergefds = {'nodes', 'faces'};

    if useVirtual
        
        prod.pivottbl = bcneumnodefacevectbl;
        
        N = bcneumnodefacetbl.num;
        [vec, i] = ind2sub([d_num, N], (1 : bcneumnodefacevectbl.num)');

        prod.dispind1 = i;
        prod.dispind2 = (1 : bcneumnodefacevectbl.num)';
        prod.dispind3 = (1 : bcneumnodefacevectbl.num)';
        
        prod.issetup = true;
        
    else
        
        prod = prod.setup();
        
    end

    extforce = extforce + prod.eval(apneum, facetNormals);

    % the format of extforce expected in assembly is in nodefacevectbl
    map = TensorMap();
    map.fromTbl  = bcneumnodefacevectbl;
    map.toTbl    = nodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual

        map.pivottbl = bcneumnodefacevectbl;

        N = bcneumnodefacetbl.num; 
        [vec, i] = ind2sub([d_num, N], (1 : bcneumnodefacevectbl.num)');

        map.dispind1 = (1 : bcneumnodefacevectbl.num)';
        map.dispind2 = sub2ind([d_num, nodefacetbl.num], vec, nodeface_from_bcneumnodeface(i));
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end

    extforce = map.eval(extforce);

    % Compute body force
    force = NaN(Nc, Nd);
    for idim = 1 : Nd
        force(:, idim) = force_fun{idim}(cc{:});
    end
    force = bsxfun(@times, G.cells.volumes, force);
    % Here, we assume we know the structure of cellvectbl;
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
    [bcnodefacetbl, indstruct] = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcface_from_bcnodeface   = indstruct{1}.inds;
    nodeface_from_bcnodeface = indstruct{2}.inds;
    bcnodefacevectbl = crossIndexArray(bcnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

    map = TensorMap();
    map.fromTbl  = nodefacevectbl;
    map.toTbl    = bcnodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual

        map.pivottbl = bcnodefacevectbl;

        N = bcnodefacetbl.num; 
        [vec, i] = ind2sub([d_num, N], (1 : bcnodefacevectbl.num)');

        map.dispind1 = sub2ind([d_num, nodefacetbl.num], vec, nodeface_from_bcnodeface(i));
        map.dispind2 = (1 : bcnodefacevectbl.num)';
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end

    bcnodefacecents = map.eval(nodefacecents);

    % Here, we assume a given structure of bcnodefacevectbl:
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
    cellvec12tbl = tbls.cellvec12tbl;
    vec12tbl = tbls.vec12tbl;

    map = TensorMap();
    map.fromTbl  = vec12tbl;
    map.toTbl    = cellvec12tbl;
    map.mergefds = {'vec1', 'vec2'};

    if useVirtual

        map.pivottbl = cellvec12tbl;

        N = tbls.celltbl.num; 
        [vec, i] = ind2sub([vec12tbl.num, N], (1 : cellvec12tbl.num)');

        map.dispind1 = vec;
        map.dispind2 = (1 : cellvec12tbl.num)';
        
        map.issetup = true;
        
    else
        
        map = map.setup();
        
    end

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
        assembly = blockAssembleBiot(G, props, drivingforces, eta, tbls, mappings, options{:}, ...
                                     'useVirtual', useVirtual);
    else
        assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, ...
                                'bcetazero', opt.bcetazero                  , ...
                                'useVirtual', useVirtual);
    end

    clear prop loadstruct

    B   = assembly.B;
    rhs = assembly.rhs;
    sol = B\rhs;

    % Displacement values at cell centers.
    cellvectbl = tbls.cellvectbl;
    ncc = cellvectbl.num;

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
