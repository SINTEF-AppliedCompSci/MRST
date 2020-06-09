function operators = setupBiotOperators(model)
    
    G = model.G;
    mech  = model.mech;
    fluid = model.fluid;
    rock  = model.rock;

    %% setup pore volume
    pv = poreVolume(G, rock);

    %% setup mechanic input
    mechprops = mech.prop;
    loadstruct = mech.loadstruct;
    
    %% setup flow input
    perm = model.rock.perm;
    
    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', false);

    celltbl = tbls.celltbl;
    colrowtbl = tbls.colrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;
    
    prod = TensorProd();
    prod.tbl1 = colrowtbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellcolrowtbl;
    prod = prod.setup();

    K = prod.eval([1; 0; 0; 1], perm);
    src = []; % no source at this stage
    
    bcstruct = fluid.bcstruct;
    
    eta = 0;
    bcetazero = false;

    %% run assembly 
    
    fluidprops.K = K;
    fluidforces.bcstruct = bcstruct;
    fluidforces.src = []; % no explicit sources here
    
    coupprops.rho = 0; % the accumulation term is set within the equation
    coupprops.alpha = rock.alpha;
    props = struct('mechprops' , mechprops , ...
                   'fluidprops', fluidprops, ...
                   'coupprops' , coupprops);
    
    drivingforces = struct('mechanics', loadstruct, ...
                           'fluid'    , fluidforces);
    
    [tbls, mappings] = setupStandardTables(G);
    
    assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, 'adoperators', true);
    
    operators = assembly.adoperators;
    operators.pv = pv;
   
end

