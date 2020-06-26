function operators = setupBiotAdOperators(model)
    
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
    
    if isfield(fluid, 'bcstruct')
        bcstruct = fluid.bcstruct;
    else
        bcstruct.bcdirichlet = [];
        bcstruct.bcneumann = [];
        bcstruct.bcneumann = [];
    end
        
    eta = model.eta;
    bcetazero = model.bcetazero;

    %% run assembly 
    
    fluidprops.K = K;
    fluidforces.bcstruct = bcstruct;
    fluidforces.src = []; % no explicit sources here (for the moment)
    
    coupprops.rho = 0; % the accumulation term is set within the equation
    coupprops.alpha = rock.alpha;
    props = struct('mechprops' , mechprops , ...
                   'fluidprops', fluidprops, ...
                   'coupprops' , coupprops);
    
    drivingforces = struct('mechanics', loadstruct, ...
                           'fluid'    , fluidforces);
    
    [tbls, mappings] = setupStandardTables(G);
    
    assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, 'addAdOperators', true);
    adoperators = assembly.adoperators;
    
    operators = setupOperatorsTPFA(G, rock);
    
    fullmomentop = adoperators.momentop;
    fulldivuop = adoperators.divuop;
    fullfacenodedispop = adoperators.facenodedispop;

    extforce = loadstruct.extforce;
    momentop = @(u, p, lm) fullmomentop(u, p, lm, extforce);
    divuop   = @(u, p, lm) fulldivuop(u, p, lm, extforce);
    facenodedispop   = @(u, p, lm) fullfacenodedispop(u, p, lm, extforce);
    
    operators.momentop        = momentop;
    operators.divuop          = divuop;
    operators.facenodedispop  = facenodedispop;
    operators.fluxop          = adoperators.fluxop;
    operators.stressop        = adoperators.stressop;
    operators.mechDirichletop = adoperators.mechDirichletop;
    
end

