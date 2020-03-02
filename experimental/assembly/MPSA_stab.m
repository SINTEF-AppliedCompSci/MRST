function out = MPSA_stab(gem, constit, varargin)
% Discretize the elliptic vector equation with given permeability field.
%
% The discretization applies a
% MPSA-0 type approximation with stabilization, as described in Nordbotten (2015/2016).
%
% Discretization should work (more or less) in both 2D and 3D.
%
% Input options:
% LSQdiv: {1} Guarantees machine precision for divD term; 0 = avoids call to lsqr (much faster)
% minWeight: {'orig'} Square root of harmonic weights in local problems; 'paper' as in m/s. 
%
% Output contains following matrices: 
% out.A - discretization of mechanics
% out.divD - divergence of displacement. Transpose can be used as gradient of pressure? 
% out.gradP - gradient of pressure
% out.stress - local stress expressions
% out.stabDelta - pressure stabilization term
% out.rhsA - RHS correction for Dirichlet boundary condition 
% out.rhsdivD - RHS correction for Dirichlet boundary condition
% out.rhsF - RHS correction for Dirichlet boundary condition
%

opt=struct('LSQdiv',1, ...
    'minWeight', 'orig');

opt=merge_options(opt,varargin{:});

LSQdiv = opt.LSQdiv; minWeight = opt.minWeight; 


Nd = gem.Nd;
Nc = gem.Nc;
Nf = gem.Nf;
Ncorn = gem.Ncorn;
warning('off', 'MATLAB:lsqr:tooSmallTolerance');

% Discretization

% Variable numbering: variable(gradient, dimension), (subcell, )cell
%
% Note: gradient is in units displacement, e.g. scaled by loc_dx
%
% System variables: (sg, sd, d)
% Full system:
% Stress match  : [ fl   0   0]  [sg] = [0]
% Edge displac. : [epg epp   0]. [sp] = [0]
% Center displ. : [cpg cpp dcp]  [ p] = [0]
% Conservation  : [con   0   0]       = [r]
%
% MPFA system: Under hypothesis/observation that upper blocks are locally invertible...
%
% [sg]   [ fl   0]-1   [     0]
% [sp] = [epg epp]   * [     0]
%        [cpg cpp]     [-dcp*p]
%
% Then the cell displacement can be found by the (relatively) small, sparse, and locally computable system matrix:
%
%             [ fl   0]-1   [   0]
% [con   0] * [epg epp]   * [   0] * [p] = [r]
%             [cpg cpp]     [-dcp]
%

% Pre-processing traction expression on faces.
temp_s2il =zeros(Nd);
for iter4 = 1:Nd
    for iter5 = 1:Nd  % 4-d tensor product
        temp_s2il(iter4,iter5) = sub2ind_loc([Nd,Nd],iter5, iter4);
    end
end

temp_val1_all = zeros(gem.Nf,Nd^2,Nd);
temp_val2_all = zeros(gem.Nf,Nd^2,Nd);
for iter3 = 1:gem.Nf
    lcell = min(gem.face_cells{iter3});
    hcell = max(gem.face_cells{iter3});
    loc_normal = gem.face_norm{iter3};
    for iter4 = 1:Nd; % Component by component....
        
        % Stress
        temp_val1 = 0;
        temp_val2 = 0;
        for iter5 = 1:Nd  % 4-d tensor product
            %temp_s2il = sub2ind_loc([Nd,Nd],iter5, iter4);
            temp_val1 = temp_val1 + constit{lcell}(temp_s2il(iter4,iter5),:)*loc_normal(iter5);
            temp_val2 = temp_val2 + constit{hcell}(temp_s2il(iter4,iter5),:)*loc_normal(iter5);
        end
        temp_val1_all(iter3,:,iter4) = temp_val1;
        temp_val2_all(iter3,:,iter4) = temp_val2;
    end
end

A = spalloc(Nc*Nd, Nc*Nd, 2 * Nc*Nd^2*3^(Nd-1));
F = spalloc(Nf*Nd, Nc*Nd, 2 * Nf*Nd*2*3^(Nd-1));
divD = spalloc(Nc, Nc*Nd, 2 * Nc*Nd*3^Nd);
gradP = spalloc(Nc*Nd, Nc, 2 * Nc*Nd*3^Nd);
stab = spalloc(Nc, Nc, 2 * Nc*Nd*3^Nd);
spLambdaIJK = zeros(Nc^2*Nd, 3); 

temp_counter_Lambda = 0; 

% Needed for non-zero Dirichlet boundary conditions.   Neglected impact of "stabilization-like" term!!!
DirichletRHS = struct();
DirichletRHS.A = spalloc(Nc*Nd, Nc*Nd, 2 * Nc*Nd^2*3^(Nd-1));
DirichletRHS.F = spalloc(Nf*Nd, Nc*Nd, 2 * Nf*Nd*2*3^(Nd-1));
DirichletRHS.divD = spalloc(Nc, Nc*Nd, 2 * Nc*Nd*3^Nd);


for iter1 = 1:Ncorn
    int_cells = gem.corner_cells{iter1};
    int_faces = gem.corner_faces{iter1};
    intmaxconstit = 0; 
    corner_count = 0;
    gco1 = gem.corner{iter1};
    stab_rhs = zeros(Nd*numel(int_faces), numel(int_cells));
    
    % Introduce characteristic face length:
    loc_dx = 0;
    
    for iter2 = int_cells'
        corner_count = corner_count + numel(gem.cell_corners{iter2});
        intmaxconstit = max(intmaxconstit, max(abs(eig(constit{iter2}))));
        loc_dx = max(loc_dx,norm(gco1-gem.center{iter2}));
    end
    
    % Something still needs to be understood here for 2D/3D
    %if ((corner_count/numel(int_cells)) > (Nd+1)) ...
    %if ((corner_count/numel(int_cells)) > (Nd+1) && numel(int_cells)>(Nd+1)) ...
    if ((corner_count/numel(int_cells)) > (Nd+1) && numel(int_cells)>(Nd+0)) ...
            || (sum(strcmp({gem.isbnd{int_faces}}, 'Neumann'))==2)
        % Standard stabilized O method - not sure why the second check is needed...
        Openalty = 1;
    else
        % Simplified integration on simplex grids.
        Openalty = 0;
    end
    Openalty=1;
    % Only perform one calculation per node
    runO = 1;
    runstopO = 0;
    
    
    
    for iter2 = int_faces'
        temp_Hookel = zeros(Nd, length(int_cells)*Nd^2);
        temp_Hookeh = zeros(Nd, length(int_cells)*Nd^2);  % Not optimal, but retains old code sturcture...
        temp_divDl = zeros(1, length(int_cells)*Nd^2);
        temp_divDh = zeros(1, length(int_cells)*Nd^2);
        
        lcell = min(gem.face_cells{iter2});
        hcell = max(gem.face_cells{iter2});
        for iter4 = 1:Nd; % Component by component....
            % Stress
            temp_val1 = temp_val1_all(iter2,:,iter4);
            temp_val2 = temp_val1_all(iter2,:,iter4);
            
            % Gradient to stress
            temp_Hookel(iter4, Nd^2*(find(lcell==int_cells)-1) + (1:Nd^2)) = temp_val1*gem.area{iter2}/numel(gem.face_corners{iter2})/loc_dx;
            temp_Hookeh(iter4, Nd^2*(find(hcell==int_cells)-1) + (1:Nd^2)) = temp_val2*gem.area{iter2}/numel(gem.face_corners{iter2})/loc_dx;
        end
        
        % Gradient to area weighted divergence
        temp_divDl(Nd^2*(find(lcell==int_cells)-1) + (1:Nd^2)) = reshape(eye(Nd), 1, Nd^2) * gem.volume{lcell}/(numel(gem.cell_faces{lcell})*numel(gem.face_corners{iter2}))/loc_dx;
        temp_divDh(Nd^2*(find(hcell==int_cells)-1) + (1:Nd^2)) = reshape(eye(Nd), 1, Nd^2) * gem.volume{hcell}/(numel(gem.cell_faces{hcell})*numel(gem.face_corners{iter2}))/loc_dx;
        
        
        if runO  % Set up and solve local problems once
            runstopO = 1;
            
            % Setup allow for "full" continuity - we will the solve lsqr for displacement
            % Matrix for traction matching
            stress_gr = zeros(length(int_faces)*Nd,length(int_cells)*Nd^2);
            flux_po = zeros(length(int_faces)*Nd,length(int_cells)*Nd);  % Remains zero
            
            % Matrices for face displacement matching at 2^(Nd-1) quadrature point
            fpot_sub_gr = zeros(length(int_faces)*Nd^2,length(int_cells)*Nd^2);
            fpot_sub_po = zeros(length(int_faces)*Nd*2^(Nd-1),length(int_cells)*Nd);
            
            % Matrices for node displacement matching
            %npot_sub_gr = zeros(length(int_cells)*Nd,length(int_cells)*Nd^2);  % Remains zero
            %npot_sub_po = eye(length(int_cells)*Nd);
            
            % RHS structures
            Dirichlet_edges = 0;
            Dirichlet_data = fpot_sub_po;
            Dirichlet_edgenum = [];
            
            eqcount_face_fl = 0;
            eqcount_face_po = 0;
            
            
            %Stress and weak displacement conditions
            for iter3 = int_faces'
                
                lcell = min(gem.face_cells{iter3});
                hcell = max(gem.face_cells{iter3});
                gfce3 = gem.face_center{iter3};
                gfn3 = gem.face_norm{iter3};
                gclc = gem.center{lcell};
                gchc = gem.center{hcell};
                dg3l = gfce3 - gclc;
                dg3h = gfce3 - gchc;
                dg3o1 = gfce3-gco1;
                if Nd == 3
                    lcdg3o1n3 = loc_cross(dg3o1, gfn3);
                end
                %loc_normal = gfn3;
                
                % Gauss points
                if Openalty && strcmp(gem.isbnd{iter3}, 'no')
                    if Nd == 2
                        eta_set = [.5-sqrt(1/3)/2, .5+sqrt(1/3)/2];
                    elseif Nd == 3
                        eta_set = [0.5-sqrt(1/3)/2, 0.5+sqrt(1/3)/2, 0.5, 0.5;
                            0, 0, -sqrt(1/3)/2, +sqrt(1/3)/2]; % Quadrature points not completely accurate in 3D (should not matter much).
                    end
                else
                    eta_set = 1/3 * ones(1,Nd-1);
                end
                
                % Harmonic weighting of displacement continuity for lsqr problem
                mean_exp = -1;
                loc_weight = ((1/2*( max(abs(eig(constit{lcell})))^mean_exp + max(abs(eig(constit{hcell})))^mean_exp ) )^(1/mean_exp)) / intmaxconstit;
                % Above as in paper. Original code uses:
                if strcmp(minWeight, 'orig')
                    loc_weight = sqrt(loc_weight);
                end
                
                ga3gfci3 = gem.area{iter3}/numel(gem.face_corners{iter3});
                for iter4 = 1:Nd; % Component by component....
                    
                    % Stress
                    temp_val1 = temp_val1_all(iter3,:,iter4);
                    temp_val2 = temp_val2_all(iter3,:,iter4);
                    
                    
                    %Rescale for proper conditioning of fine-scale problems
                    temp_val1 = temp_val1/intmaxconstit / loc_dx;
                    temp_val2 = temp_val2/intmaxconstit / loc_dx;
                    
                    
                    if ~strcmp(gem.isbnd{iter3}, 'Dirichlet')
                        eqcount_face_fl=eqcount_face_fl+1;
                        
                        stress_gr(eqcount_face_fl, Nd^2*(find(lcell==int_cells)-1) + (1:Nd^2)) = temp_val1*ga3gfci3;
                        stab_rhs(eqcount_face_fl, find(lcell==int_cells)) = -gfn3(iter4)*ga3gfci3 /intmaxconstit / loc_dx;
                        
                        if ~strcmp(gem.isbnd{iter3}, 'Neumann')  
                            stress_gr(eqcount_face_fl, Nd^2*(find(hcell==int_cells)-1) + (1:Nd^2)) = - temp_val2*ga3gfci3;
                            stab_rhs(eqcount_face_fl, find(hcell==int_cells)) = gfn3(iter4)*ga3gfci3 /intmaxconstit / loc_dx;
                        end
                    end
                    
                    % Continuity at Gauss continuity points for displacement
                    if ~strcmp(gem.isbnd{iter3}, 'Neumann')
                        for iter6 = 1:length(eta_set(1,:));
                            eqcount_face_po=eqcount_face_po+1;
                            
                            eta = eta_set(:,iter6);
                            
                            loc_x = dg3l - eta(1)*dg3o1;
                            neigh_loc_x = dg3h - eta(1)*dg3o1;
                            
                            if length(eta)==2
                                loc_x = loc_x + eta(2)*lcdg3o1n3;
                                neigh_loc_x = neigh_loc_x + eta(2)*lcdg3o1n3;
                            end
                            
                            
                            fpot_sub_gr(eqcount_face_po, Nd^2*(find(hcell==int_cells)-1) + Nd*(iter4-1) + (1:Nd)) = -loc_weight*neigh_loc_x/loc_dx;
                            fpot_sub_po(eqcount_face_po, Nd*(find(hcell==int_cells)-1)+iter4) = -loc_weight*1;
                            
                            if strcmp(gem.isbnd{iter3}, 'Dirichlet')
                                Dirichlet_edges = Dirichlet_edges + 1;
                                Dirichlet_data(eqcount_face_po, Dirichlet_edges) = 1;
                                Dirichlet_edgenum = [Dirichlet_edgenum; iter3];
                            elseif strcmp(gem.isbnd{iter3}, 'no')
                                fpot_sub_gr(eqcount_face_po, Nd^2*(find(lcell==int_cells)-1) + Nd*(iter4-1) + (1:Nd)) = loc_weight*loc_x/loc_dx;
                                fpot_sub_po(eqcount_face_po, Nd*(find(lcell==int_cells)-1)+iter4) = loc_weight*1;
                            end
                        end
                    end
                end
            end
            
            % Only solve local problems once [could be structured better]
            
            % Remove reduntant lines
            stress_gr = stress_gr(1:eqcount_face_fl,:);
            flux_po = flux_po(1:eqcount_face_fl,:);
            fpot_sub_gr = fpot_sub_gr(1:eqcount_face_po,:);
            fpot_sub_po = fpot_sub_po(1:eqcount_face_po,:);
            
            if Dirichlet_edges > 0
                Dirichlet_data = Dirichlet_data(1:eqcount_face_po,1:Dirichlet_edges);
            end
            
            disc_gr = zeros(length(int_cells)*Nd,length(int_cells)*Nd^2);
            disc_po = eye(length(int_cells)*Nd);
            
            %Assemble local system
            local_matrices = [  stress_gr  ,          flux_po; ...
                fpot_sub_gr,        fpot_sub_po; ...
                disc_gr,            disc_po            ];
            
            %Nd RHS per node center.
            npot_disc = speye(length(int_cells)*Nd);
            rhs = [zeros(eqcount_face_fl+eqcount_face_po,length(int_cells)*Nd); npot_disc];
            
            % Additional 1 RHS per node center (only for use with Biot!)
            RHS = {rhs, [stab_rhs * loc_dx; zeros(eqcount_face_po+length(int_cells)*Nd, length(int_cells))]};
            
            for iter9 = 1:2;
                if (Openalty~=1)
                    local_solve{iter9} = local_matrices \ RHS{iter9};
                else
                    [~, slm2]  = size(local_matrices);
                    
                    % Minimize displacement jumps
                    minvars = (eqcount_face_fl+1):(eqcount_face_fl+eqcount_face_po);
                    % Enforce continuity of strictly
                    constraints = setdiff(1:length(local_matrices), minvars);
                    
                    local_square = local_matrices(minvars,:)' * local_matrices(minvars,:);
                    rhs_square = local_matrices(minvars,:)' * RHS{iter9}(minvars,:);
                    
                    local_constraints = local_matrices(constraints, :);
                    local_matrices2 = [local_square, local_constraints'; ...
                        local_constraints, zeros(length(constraints))];
                    rhs2 = [rhs_square; RHS{iter9}(constraints,:)];
                    
                    lc_1 = local_constraints(1:eqcount_face_fl,:);
                    [s1, s2, s3] = svd(lc_1');
                    m1 = find(abs(diag(s2))>10*eps*max(diag(s2)));
                    if length(m1)==eqcount_face_fl
                        % Direct solve.
                        local_solve2 = local_matrices2 \ rhs2;
                    elseif numel(int_cells) < Nd % Not tested for Nd == 3
                        s4 = s2; 
                        s4((length(m1)+1):Nd^2,(length(m1)+1):Nd^2) = eye(Nd^2-length(m1));
                        lc_2 = s1*s4*s3'; 
                        local_matrices3 = local_matrices;
                        local_matrices3(1:eqcount_face_fl,:) = lc_2'; 
                        local_solve2 = local_matrices3 \ RHS{iter9};
                    else
                        % Eliminate redundant stress conditions (1(6) for orthogonal cart. grids in 2D (3D), variable in general)
                        lc_2 = s3(:,m1)'*lc_1;
                        local_constraints2 = [lc_2; local_constraints((eqcount_face_fl+1):end,:)];
                        local_matrices3 = [local_square, local_constraints2'; ...
                            local_constraints2, zeros(size(local_constraints2, 1))];
                        rhs3 = [rhs_square; s3(:,m1)'*RHS{iter9}(constraints(1:eqcount_face_fl),:); RHS{iter9}(constraints((eqcount_face_fl+1):end),:)];
                        local_solve2 = local_matrices3 \ rhs3;
                        
                        if LSQdiv
                            % Correct solution so divD is consistent for constants to machine precision
                            for iter3 = 1:size(RHS{iter9},2)
                                [lloc_inv, ~] = lsqr(local_matrices3, rhs3(:,iter3), 10^-100, 1, [], [], local_solve2(:,iter3));
                                %[lloc_inv, FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(local_matrices3, rhs3(:,iter3), 10^-100, 1, [], [], local_solve2(:,iter3));
                                local_solve2(:,iter3) = lloc_inv;
                            end
                        end
                    end
                    local_solve{iter9} = local_solve2(1:slm2, :);
                end
                
                if runstopO == 1
                    runO = 0;
                end
            end
        end
        
        
        % Extract displacement and calculate divergence
        for iter9 = 1:2
            f_sten = temp_Hookel*local_solve{iter9}(1:(length(int_cells)*Nd^2), :);
            f_stenl{iter9} = f_sten.*(abs(f_sten)>10*eps*max(max(f_sten)));
            f_sten = temp_Hookeh*local_solve{iter9}(1:(length(int_cells)*Nd^2), :);
            f_stenh{iter9} = f_sten.*(abs(f_sten)>10*eps*max(max(f_sten)));
            
            divD_sten = temp_divDl*local_solve{iter9}(1:(length(int_cells)*Nd^2), :);
            divD_stenl{iter9} = divD_sten.*(abs(divD_sten)>10*eps*max(max(divD_sten)));
            divD_sten = temp_divDh*local_solve{iter9}(1:(length(int_cells)*Nd^2), :);
            divD_stenh{iter9} = divD_sten.*(abs(divD_sten)>10*eps*max(max(divD_sten)));
        end
        
        if Dirichlet_edges > 0
            Dirichlet_rhs = [zeros(eqcount_face_fl, Dirichlet_edges); Dirichlet_data; zeros(length(int_cells)*Nd,Dirichlet_edges)];
            Dirichlet_local_solve = local_matrices \ Dirichlet_rhs;
            
            % FV type boundary treatment
            Dirichlet_f_sten = temp_Hookel*Dirichlet_local_solve(1:(length(int_cells)*Nd^2), :);
            Dirichlet_f_sten = Dirichlet_f_sten.*(abs(Dirichlet_f_sten) > (10*eps*max(max(Dirichlet_f_sten))));
            
            Dirichlet_divD_stenl = temp_divDl*Dirichlet_local_solve(1:(length(int_cells)*Nd^2), :);
            Dirichlet_divD_stenl = Dirichlet_divD_stenl.*(abs(Dirichlet_divD_stenl)>10*eps*max(max(Dirichlet_divD_stenl)));
            Dirichlet_divD_stenh = temp_divDh*Dirichlet_local_solve(1:(length(int_cells)*Nd^2), :);
            Dirichlet_divD_stenh = Dirichlet_divD_stenh.*(abs(Dirichlet_divD_stenh)>10*eps*max(max(Dirichlet_divD_stenh)));
        end
        
        
        if iter1 == 2
            iter1;
        end
        
        loc_cells = gem.face_cells{iter2};
        temp_vec = zeros(length(int_cells)*Nd,1);
        for iter6 = 1:Nd
            temp_vec(iter6:Nd:end) = Nd*(int_cells-1)+iter6;
        end
        
        % Add incoming fluxes
        A(Nd*(min(loc_cells)-1)+(1:Nd),temp_vec) = A(Nd*(min(loc_cells)-1)+(1:Nd),temp_vec) + f_stenl{1};
        % Add contribution to divergence of displacement
        divD(min(loc_cells),temp_vec) = divD(min(loc_cells),temp_vec) + divD_stenl{1};
        % Add contribution to gradient of pressure
        gradP(Nd*(min(loc_cells)-1)+(1:Nd),int_cells) = gradP(Nd*(min(loc_cells)-1)+(1:Nd),int_cells) + f_stenl{2};
        % Add contribution to stabilization
        stab(min(loc_cells),int_cells) = stab(min(loc_cells),int_cells) - divD_stenl{2};
        
        if numel(loc_cells) == 2
            A(Nd*(max(loc_cells)-1)+(1:Nd),temp_vec) = A(Nd*(max(loc_cells)-1)+(1:Nd),temp_vec) - f_stenl{1};  %Reuse f_stenl since forces are equal and opposite
            divD(max(loc_cells),temp_vec) = divD(max(loc_cells),temp_vec) + divD_stenh{1};
            gradP(Nd*(max(loc_cells)-1)+(1:Nd),int_cells) = gradP(Nd*(max(loc_cells)-1)+(1:Nd),int_cells) - f_stenh{2}; %Note use of f_stenh since forces not equal!!!
            stab(max(loc_cells),int_cells) = stab(max(loc_cells),int_cells) - divD_stenh{2};
        end
        
        F(Nd*(iter2-1)+(1:Nd),temp_vec) = F(Nd*(iter2-1)+(1:Nd),temp_vec)+f_stenl{1};
        % Here it is possible to add the pressure influence on surface stress, similar to the stab term. May not be important, but maybe.
        
        
        %            spLambdaIJK(temp_counter_Lambda +(1:Nd),1) = (int_cells(iter4)-1)*Nd + (1:Nd);
        %            spLambdaIJK(temp_counter_Lambda +(1:Nd),2) = int_cells(iter5);
        %            spLambdaIJK(temp_counter_Lambda +(1:Nd),3) = temp_grad(:,iter5)'*temp_prod(:,(iter4-1)*Nd +(1:Nd)) * (gem.volume{iter6}/(numel(gem.cell_corners{iter6})));
        %            temp_counter_Lambda = temp_counter_Lambda + Nd;
        
        
        for iter6 = 1:(Dirichlet_edges/Nd)
            % Cell inside boundary edge
            loc_cell2 = gem.face_cells{Dirichlet_edgenum(Nd*iter6)};
            
            % Contribution to mass balance of neighbours to iter2 due to bnd.condition at iter6.
            DirichletRHS.A(Nd*(min(loc_cells)-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd)) = DirichletRHS.A(Nd*(min(loc_cells)-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd)) + Dirichlet_f_sten(:,Nd*(iter6-1)+(1:Nd));
            DirichletRHS.divD(min(loc_cells),Nd*(loc_cell2-1)+(1:Nd)) = DirichletRHS.divD(min(loc_cells),Nd*(loc_cell2-1)+(1:Nd)) + Dirichlet_divD_stenl(:,Nd*(iter6-1)+(1:Nd));
            
            % Should a term like this exist?
            % DirichletRHS.gradP = ?
            
            if numel(loc_cells) == 2
                DirichletRHS.A(Nd*(max(loc_cells)-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd)) = DirichletRHS.A(Nd*(max(loc_cells)-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd)) - Dirichlet_f_sten(:,Nd*(iter6-1)+(1:Nd));
                DirichletRHS.divD(max(loc_cells),Nd*(loc_cell2-1)+(1:Nd)) = DirichletRHS.divD(max(loc_cells),Nd*(loc_cell2-1)+(1:Nd)) + Dirichlet_divD_stenh(:,Nd*(iter6-1)+(1:Nd));
            end
            
            DirichletRHS.F(Nd*(iter2-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd)) = DirichletRHS.F(Nd*(iter2-1)+(1:Nd),Nd*(loc_cell2-1)+(1:Nd))+Dirichlet_f_sten(:,Nd*(iter6-1)+(1:Nd));
        end
        
    end
end


warning('off', 'MATLAB:lsqr:tooSmallTolerance');

out.A = A; clear A;
out.divD = divD; clear divD;
out.gradP = gradP; clear gradP;
out.stress = F; clear F;
out.stabDelta = stab; clear stab;
%out.stabLambda = sparse(spLambdaIJK(:,1),spLambdaIJK(:,2),spLambdaIJK(:,3));
out.rhsA = DirichletRHS.A;
out.rhsdivD = DirichletRHS.divD;
out.rhsF = DirichletRHS.F;
clear DirichletRHS;

end

function c = loc_cross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
    a(3)*b(1)-a(1)*b(3);
    a(1)*b(2)-a(2)*b(1)];
end