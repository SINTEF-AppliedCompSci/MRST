function assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, varargin)
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


    opt = struct('bcetazero'       , true, ...
                 'assemblyMatrices', false, ...
                 'addAdOperators'  , false, ...
                 'scalingfactor', []);
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

    if ~isempty(opt.scalingfactor)
        scalingfactor = opt.scalingfactor;
        usescaling = true;
    else
        usescaling = false;
    end
    
    mechprops  = props.mechprops;
    fluidprops = props.fluidprops;
    coupprops  = props.coupprops;
    
    % Assemble mechanic problem
    loadstruct = drivingforces.mechanics;
    options = {'assemblymatrices', true, ...
               'bcetazero'       , opt.bcetazero};
    mechassembly = assembleMPSA(G, mechprops, loadstruct, eta, tbls, mappings, options{:});
    
    % Assemble fluid problem
    fluidforces = drivingforces.fluid;
    bcstruct = fluidforces.bcstruct;
    src = fluidforces.src;
    K = fluidprops.K;
    options = {'addAdOperators', true, ...
               'bcetazero'     , opt.bcetazero};
    fluidassembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, options{:});
    
    % Assemble coupling terms (finite volume and consistent divergence operators)
    % we need to pass the number of nodes per cell to assembleCouplingTerms
    alpha = coupprops.alpha;

    cellnodetbl = tbls.cellnodetbl;
    celltbl = tbls.celltbl;
    cell_from_cellnode = mappings.cell_from_cellnode;
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = celltbl;     
    map.mergefds = {'cells'};

    map.pivottbl = cellnodetbl;
    ind = (1 : cellnodetbl.num)';
    map.dispind1 = ind;
    map.dispind2 = cell_from_cellnode(ind);
    map.issetup = true;
    
    nnodespercell = map.eval(ones(cellnodetbl.num, 1));

    coupassembly = assembleCouplingTerms(G, eta, alpha, nnodespercell, tbls, mappings);
    
    % Recover matrices from mechanic assembly
    mechmat = mechassembly.matrices;
    invA11 = mechmat.invA11;
    A11 = mechmat.A11;
    A12 = mechmat.A12;
    A21 = mechmat.A21;
    A22 = mechmat.A22;
    A15 = -mechmat.D;
    A51 = -A15';
    C1  = mechmat.C1;
    C2  = mechmat.C2;
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
    A14 = -A14'; % We use the gradient which is the transpose of minus div
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
    %       |      B42       B44 |
    
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
    
    fac = 1e-9;
    if usescaling
        B1 = scalingfactor*B1;
        B2 = scalingfactor*B2;
    end
    
    B = [B1; B2; B3; B4];
    
    % Assembly of right hand side
    f = fullrhs; % shortcut
    redrhs = cell(4, 1);
    redrhs{1} = f{2} - A21invA11*f{1};
    redrhs{2} = f{4} - A41invA11*f{1} - A43invA33*f{3};
    redrhs{3} = f{5} - A51invA11*f{1};
    redrhs{4} = f{6} - A63invA33*f{3};
    
    if usescaling
        redrhs{1} = scalingfactor*redrhs{1};
        redrhs{2} = scalingfactor*redrhs{2};
    end
    rhs = vertcat(redrhs{:});

    assembly = struct('B'  , B, ...
                      'rhs', rhs);
   
    if opt.assemblyMatrices
        fullsystem.A = A;
        fullsystem.rhs = vertcat(fullrhs{:});
        assembly.fullsystem = fullsystem;
    end    
    
    if opt.addAdOperators

        fluxop = fluidassembly.adoperators.fluxop;
        
        % Setup face node displacement operator
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
