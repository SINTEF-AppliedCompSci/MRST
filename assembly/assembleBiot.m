function assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, varargin)
    
    opt = struct('assemblyMatrices', false, ...
                 'adoperators'     , false);
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
        
    mechprops  = props.mechprops;
    fluidprops = props.fluidprops;
    coupprops  = props.coupprops;
    
    % Assemble mechanic problem
    loadstruct = drivingforces.mechanics;
    mechassembly = assembleMPSA(G, mechprops, loadstruct, eta, tbls, mappings, 'assemblymatrices', true);
    
    % Assemble fluid problem
    fluidforces = drivingforces.fluid;
    bcstruct = fluidforces.bcstruct;
    src = fluidforces.src;
    K = fluidprops.K;
    fluidassembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings);

    % Assemble coupling terms (finite volume and consistent divergence operators)
    alpha = coupprops.alpha;
    coupassembly = assembleCouplingTerms(G, eta, alpha, tbls, mappings);
    
    % Recover matrices from mechanic assembly
    mechmat = mechassembly.matrices;
    invA11 = mechmat.invA11;
    A11 = mechmat.A11;
    A12 = mechmat.A12;
    A21 = mechmat.A21;
    A22 = mechmat.A22;
    A15 = -mechmat.D;
    A51 = -A15';
    mechrhs = mechmat.fullrhs;
    
    % Recover matrices from fluid assembly
    fluidmat = fluidassembly.matrices;
    invA33 = fluidmat.invA11;
    A33 = fluidmat.A11;
    A34 = fluidmat.A12;
    A43 = fluidmat.A21;
    A44 = fluidmat.A22;
    A36 = -fluidmat.D;
    A63 = -A36';
    fluidrhs = fluidmat.fullrhs;
    
    % Recover the coupling terms
    A14 = coupassembly.divfv;
    A14 = -A14'; % we use the gradient which is the transpose of minus div
    A41 = coupassembly.divconsnf;
    A42 = coupassembly.divconsc;
    
    % We add the diagonal term for the mass conservation equation
    rho = coupprops.rho;
    rho = G.cells.volumes.*rho;
    
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
    fullrhs{1} = mechrhs{1};
    fullrhs{2} = mechrhs{2};
    fullrhs{3} = fluidrhs{1};
    fullrhs{4} = fluidrhs{2};
    fullrhs{5} = mechrhs{3};
    fullrhs{6} = fluidrhs{3};
    
    %    
    %        | A11    A12    0      A14    A15    0    |
    %        | A21    A22    0      0      0      0    |
    %  A =   | 0      0      A33    A34    0      A36  |
    %        | A41    A42    A43    A44    0      0    |
    %        | A51    0      0      0      0      0    |
    %        | 0      0      A63    0      0      0    | 
    %
    
    n1 = tbls.nodefacecoltbl.num;
    n2 = tbls.cellcoltbl.num;
    n3 = tbls.nodefacetbl.num;
    n4 = tbls.celltbl.num;
    n5 = size(mechrhs{3}, 1);
    n6 = size(fluidrhs{3}, 1);

    if opt.assemblyMatrices
        A1 = [A11, A12, zeros(n1, n3), A14, A15, zeros(n1, n6)];
        A2 = [A21, A22, zeros(n2, n3 + n4 + n5 + n6)];    
        A3 = [zeros(n3, n1 + n2), A33, A34, zeros(n3, n5), A36];
        A4 = [A41, A42, A43, A44, zeros(n4, n5 + n6)];
        A5 = [A51, zeros(n5, n2 + n3 + n4 + n5 + n6)];
        A6 = [zeros(n6, n1 + n2), A63, zeros(n6, n4 + n5 + n6)];
        
        A = [A1; A2; A3; A4; A5; A6];
    end
    
    %       | B11  B12  B13      |
    %  B =  | B21  B22  B23  B24 |
    %       | B31  B32  B33      |
    %       |      B32       B44 |
    
    % 1. row : Momentum equation
    A21invA11 = A21*invA11;
    B11 = -A21invA11*A12 +  A22;
    B12 = -A21invA11*A14;
    B13 = -A21invA11*A15;

    % 2. row : Fluid mass conservation
    A41invA11 = A41*invA11;
    A43invA33 = A43*invA33;
    B21 = -A41invA11*A12 + A42;
    B22 = -A41invA11*A14 - A43invA33*A34 + A44;
    B23 = -A41invA11*A15;
    B24 = -A43invA33*A36;

    % 3. row : Mechanic BC
    A51invA11 = A51*invA11;
    B31 = -A51invA11*A12;
    B32 = -A51invA11*A14;
    B33 = -A51invA11*A15;

    % 4. row : Fluid BC
    A63invA33 = A63*invA33;
    B42 = -A63invA33*A34;
    B44 = -A63invA33*A36;

    % Assemble B matrix
    n1 = tbls.cellcoltbl.num;
    n2 = tbls.celltbl.num;
    n3 = size(mechrhs{3}, 1);
    n4 = size(fluidrhs{3}, 1);
    B1 = [B11, B12, B13, zeros(n1, n4)];
    B2 = [B21, B22, B23, B24];
    B3 = [B31, B32, B33, zeros(n3, n4)];
    B4 = [zeros(n4, n1), B42, zeros(n4, n3), B44];
    
    B = [B1; B2; B3; B4];
    
    % Assembly of right hand side
    f = fullrhs; % shortcut
    rhs = cell(4, 1);
    rhs{1} = f{2} - A21invA11*f{1};
    rhs{2} = f{4} - A41invA11*f{1} - A43invA33*f{3};
    rhs{3} = f{5} - A51invA11*f{1};
    rhs{4} = f{6} - A63invA33*f{3};
    
    rhs = vertcat(rhs{:});

    assembly = struct('B'  , B, ...
                      'rhs', rhs);
   
    if opt.assemblyMatrices
        fullsystem.A = A;
        fullsystem.rhs = cat(fullrhs{:});
        assembly.fullsystem = fullsystem;
    end    
    
end

