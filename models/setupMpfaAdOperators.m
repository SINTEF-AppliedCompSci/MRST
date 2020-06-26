function operators = setupMpfaAdOperators(model)    
    
    %% We setup the usual operators and the mpfa flux operator
    %
    % we set up the mappings
    %
    %   F1 : nodefacetbl -> intfacetbl 
    %   F2 : celltbl -> intfacetbl
    %
    %  such that the flux u in intfacetbl (interior faces) is given by
    %
    %  u  = [F1  F2] * [ pnf (pressure at nodefacetbl);
    %                    pc  (pressure at celltbl)];
    %
    %  Then, we proceed with the reduction to remove dependency in pnf (pressure at nodefacetbl)
    %
    
    G         = model.G;
    rock      = model.rock;
    eta       = model.eta;
    bcetazero = model.bcetazero;
    
    perm = rock.perm;
    
    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');
    
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', false);

    celltbl = tbls.celltbl;
    colrowtbl = tbls.colrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;%
    
    prod = TensorProd();
    prod.tbl1 = colrowtbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellcolrowtbl;
    prod = prod.setup();

    K = prod.eval([1; 0; 0; 1], perm);
    src = []; % no source at this stage
    bcstruct.bcdirichlet = []; % no Dirichlet
    bcstruct.bcneumann = []; % zero neumann bc
    
    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'onlyAssemblyMatrices', true);
    mpfaKgrad = setupMpfaFlux(G, assembly, tbls);
    
    operators = setupOperatorsTPFA(G, rock);
    operators.mpfaKgrad = mpfaKgrad;
    
end

