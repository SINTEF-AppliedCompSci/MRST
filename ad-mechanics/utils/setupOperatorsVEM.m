function operators = setupOperatorsVEM(G, C, el_bc, load, alpha_scaling, S)
    
    [~, extra] = VEM_linElast(G                    , ...
                              C                   , ...
                              el_bc               , ...
                              load                , ...
                              'alpha_scaling', alpha_scaling , ...
                              'S', S                         , ...
                              'linsolve', @(A, rhs) 0 * rhs); 
    
    operators.mech = extra.disc; 
    
    vdiv    = VEM_div(G); 
    [~, op] = VEM_assemble(G, C);%, 'blocksize', model.GMech.cells.num/10);
    strain  = op.WC' * op.assemb'; 
    stress  = op.D * strain; %op.WC' * op.assemb'; 
    
    operators.extra = struct('vdiv', vdiv, 'stress', stress, 'strain', strain); 
        
    
end