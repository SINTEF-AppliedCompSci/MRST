function mpsaPaperConvergenceTest(Nd, nref)
% Nd : spatial dimension
% nref : number of refinement level
    
% Convergence test case
% title={Finite volume methods for elasticity with weak symmetry},
% author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
% journal={International Journal for Numerical Methods in Engineering},
% volume={112},
% number={8},
% pages={939--962},
% year={2017},
% publisher={Wiley Online Library}

    mrstModule add vemmech mpsaw vem


    % Number of refinement levels
    nref = 4; 

    % Test case number (see definitions below)
    testCase = 1; 

    % The constant eta (between 0 and 1) in MPSAW which defines the position of
    % continuity point

    switch testCase
      case 1
        % Cartesian grid
        gridType = 1; 
        eta = 1/3;
      case 2
        % Cartesian grid
        gridType = 1; 
        eta = 0; 
      case 3
        % Triangular grid, 90 degree angles
        gridType = 2; 
        eta = 1/3; 
      case 4
        % Equilateral triangles
        gridType = 3; 
        eta = 1/3;
      otherwise
        error('testCase not recognized');
    end

    kappa = 1;
    alpha = 1;
    % the coefficients mu and lambda are defined through the coefficients kappa
    % and alpha

    [u, f] = analyticalReferencePaper(Nd);

    for iter1 = 1 : nref
        
        disp(['Refinement ' num2str(iter1)])
        
        Nx = 2^iter1*ones(1, Nd); 
        [G, cornercelltbl] = gridForConvTest(Nx, gridType); 
        G = computeVEMGeometry(G); 
        G = computeGeometryCalc(G); 
        [tbls, mappings] = setupStandardTables(G);

        Nc = G.cells.num; 
        Nf = G.faces.num; 
        Nn = G.nodes.num; 
        Nd = G.griddim; 

        mu = ones(Nc, 1);
        ccells = cornercelltbl.get('cells');
        mu(ccells) = kappa;
        
        lambda = alpha*mu;
        
        isBoundary = any(G.faces.neighbors == 0, 2); 
        bcfaces =  find(isBoundary);
        
        
        
        xf = G.faces.centroids; 
        xc = G.cells.centroids; 
        
        clear prop
        prop.mu = mu*ones(Nc, 1); 
        prop.lambda = lambda*ones(Nc, 1); 
        
        
        clear bc
        bc{1}.extfaces = bcfaces;
        bc{1}.linform = [1; 0];
        bc{2}.extfaces = bcfaces;
        bc{2}.linform = [0; 1];
        
        clear loadstruct
        loadstruct.bc = bc;
        loadstruct.force = rhsMech;
        loadstruct.extforce = zeros(tbls.nodefacecoltbl.num, 1);
        
        assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings);
        
        B   = assembly.B  ;
        rhs = assembly.rhs;
        sol = B\rhs;

        % displacement values at cell centers.
        cellcoltbl = tbls.cellcoltbl;
        n = cellcoltbl.num;

        u = sol(1 : n);
        u = reshape(u, 2, [])';
        
        dnum = u;
        
        san1 = sum([s11(xf( :, 1), xf( :, 2)), s21(xf( :, 1), xf( :, 2))] .* G.faces.normals, 2); 
        san2 = sum([s12(xf( :, 1), xf( :, 2)), s22(xf( :, 1), xf( :, 2))] .* G.faces.normals, 2); 
        
        toc; 
        % Analytical solution
        dex = [d1(xc( :, 1), xc( :, 2)) d2(xc( :, 1), xc( :, 2))]; 
        
        % Errors in L2 and max - norm
        deL2(iter1) = sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, (dex - dnum).^2)))) / sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, ( dex).^2)))); 
        dem(iter1) = max(max(abs(dex - dnum))); 
        
        dostresscomputation = false;
        if dostresscomputation
            stress = md.stress*reshape(dnum', [], 1); 
            
            s_ex = reshape([san1, san2]', [], 1); 
            
            sem(iter1) = max(abs(s_ex - stress))/ max(abs(stress)); 
            fa = reshape(repmat(G.faces.areas, 1, Nd)', [], 1); 
            seL2(iter1) = sqrt(sum(fa.^2 .* (stress - s_ex).^2)) / sqrt(sum(fa.^2 .* s_ex.^2)); 
        end
        % make VEM solution
        
        [E, nu] = elasticModuloTransform(lambda, mu, 'lam_mu', 'E_nu'); 
        %{
        if(G.griddim == 2)
            [E, nu] = LMu2ENu_2D(lambda, mu); 
            else
                [E, nu] = LMu2ENu_3D(lambda, mu); 
                end
                %} 
                Ev = repmat(E, G.cells.num, 1); 
                nuv = repmat(nu, G.cells.num, 1); 
                C = Enu2C(Ev, nuv, G); 
                %}
                %{
                C = nan(G.cells.num, numel(constit{1})); 
                for i = 1 : G.cells.num
                    C(i, : ) = constit{i}( : ); 
                    end
                    %}
                    % set all boundary to no displacement
                    faces = find(any(G.faces.neighbors == 0, 2)); 
                    inodes = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 1); 
                    nodes = unique(G.faces.nodes(inodes)); 
                    el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', zeros(numel(nodes), G.griddim), 'faces', faces,...
                                                     'uu_face', zeros(numel(nodes), G.griddim), 'mask', true(numel(nodes), G.griddim)),...
                                   'force_bc', []); 
                    load = @(coord) - [mrhs1(coord( :, 1), coord( :, 2)), mrhs2(coord( :, 1), coord( :, 2))]; 
                    
                    tic; [uVEM, extra] = VEM_linElast(G, C, el_bc, load); toc; 
                    deVEM(iter1) = sqrt(sum(sum(bsxfun(@times, (uVEM(G.cells.nodes, : ) - dvec(G.nodes.coords(G.cells.nodes, : ))).^2, G.weights.cell_nodes), 2)))./...
                        sqrt(sum(sum(bsxfun(@times, dvec(G.nodes.coords(G.cells.nodes, : )).^2, G.weights.cell_nodes), 2))); 
                    % make CC solution with global interface
                    tic; [uu, out] = CC_linElast(G, C, el_bc, load); toc; 
                    
                    
    end

    %% Print convergence rates

    log2(deL2(1 : end - 1)./deL2(2 : end))
    log2(deVEM(1 : end - 1)./deVEM(2 : end))

    %% Error at last iteration

    fprintf('relative L2 error exact vs mpsa: %g\n', deL2(end));
    fprintf('relative L2 error exact vs vem: %g\n', deVEM(end));


    %% 
    n = 3; 
    figure
    set(gcf, 'numbertitle', 'off', 'name', 'DISPLACEMENT')

    subplot(3, 2, 1)
    title('x-mpsa')
    plotCellData(G, dnum( :, 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 2)
    title('y-mpsa')
    plotCellData(G, dnum( :, 2), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 3)
    title('x-vem')
    plotNodeData(G, uVEM( :, 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 4)
    title('y-vem')
    plotNodeData(G, uVEM( :, 2), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 5)
    title('x-exact')
    plotCellData(G, dex( :, 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 6)
    title('y-exact')
    plotCellData(G, dex( :, 2), 'edgecolor', 'none'), colorbar

    return
    %% 
    figure
    set(gcf, 'numbertitle', 'off', 'name', 'ERROR')
    uu_nn = dvec(G.nodes.coords); 
    uu_cc = dvec(G.cells.centroids); 
    val = uVEM - uu_nn; 
    subplot(2, 2, 1), 
    plotNodeData(G, val( :, 1)); colorbar
    title('x-vem')
    subplot(2, 2, 2), 
    plotNodeData(G, val( :, 2));colorbar
    title('y-vem')
    val = dnum - uu_cc;
    subplot(2, 2, 3), 
    plotCellData(G, val( :, 1));colorbar
    title('x-mpsa')
    subplot(2, 2, 4), 
    plotCellData(G, val( :, 2));colorbar
    title('x-mpsa')

    %%
    figure
    subplot(2, 1, 1), 
    plotCellData(G, mrhs1(G.cells.centroids( : , 1), G.cells.centroids( :, 2)));colorbar
    subplot(2,1,2),
    plotCellData(G,mrhs2(G.cells.centroids(:,1),G.cells.centroids(:,2)));colorbar
end


function y = Xindfunction(x)
    ix = all(x > 0.5, 2);
    y = 0*x;
    y(ix) = 1;
end
