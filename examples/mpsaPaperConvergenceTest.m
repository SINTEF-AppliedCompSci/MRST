function mpsaPaperConvergenceTest(Nd, nref, kappa, alpha)
% Nd : spatial dimension
% nref : number of refinement level
% kappa : coefficient used for the top-corner (kappa = 1 corresponds to
%         homogeneous case)
% alpha : coefficient defining lambda from mu, lambda = alpha*mu;
    
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

    % Test case number (see definitions below)
    testCase = 1; 
    doVem = false;
    
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


    [u_fun, force_fun, mu_fun] = analyticalReferencePaper(Nd, kappa);

    for iter1 = 1 : nref
        
        disp(['Refinement ' num2str(iter1)])
        
        Nx = 2^iter1*ones(1, Nd); 
        G = gridForConvTest(Nx, gridType); 
        G = computeVEMGeometry(G); 
        G = computeGeometryCalc(G); 
        [tbls, mappings] = setupStandardTables(G);
        coltbl = tbls.coltbl;
        nodefacetbl = tbls.nodefacetbl;
        nodefacecoltbl = tbls.nodefacecoltbl;
        
        Nc = G.cells.num; 
        Nf = G.faces.num; 
        Nn = G.nodes.num; 
        Nd = G.griddim; 
        
        % prepare input for analytical functions
        for idim = 1 : Nd
            cc{idim} = G.cells.centroids(:, idim);
        end
        mu = mu_fun(cc{:});
        lambda = alpha*mu;
        
        prop.mu = mu; 
        prop.lambda = lambda; 
        
        isBoundary = any(G.faces.neighbors == 0, 2); 
        bcfaces =  find(isBoundary);
        
        bcfacetbl.faces = bcfaces;
        bcfacetbl = IndexTable(bcfacetbl);
        bcnodefacetbl = crossTable(bcfacetbl, nodefacetbl, {'faces'});
        bcnodefacecoltbl = crossTable(bcnodefacetbl, coltbl, {}, 'optpureproduct', ...
                                      true);
        clear bcfacetbl
        
        [~, nodefacecents] = computeNodeFaceCentroids(G, tbls, eta);
        % Here, we assume a given structure of nodefacecoltbl:
        nodefacecents = reshape(nodefacecents, Nd, [])';
        
        map = TensorMap();
        map.fromTbl = nodefacecoltbl;
        map.toTbl = bcnodefacecoltbl;
        map.mergefds = {'nodes', 'faces', 'coldim'};
        map = map.setup();
        
        bcnodefacecents = map.eval(nodefacecents);
        % Here, we assume a given structure of bcnodefacecoltbl:
        bcnodefacecents = reshape(bcnodefacecents, Nd, [])';
        bcnum = bcnodefacetbl.num;
        
        % Prepare input for analytical functions
        for idim = 1 : Nd
            bnfc{idim} = bcnodefacecents(:, idim);
        end
        
        % Compute boundary conditions
        for idim = 1 : Nd
            linform = zeros(bcnum, Nd);
            linform(:, idim) = 1;
            linforms{idim} = linform;
            linformvals{idim} = u_fun{idim}(bnfc{:});
        end
        
        bcfaces = bcnodefacetbl.get('faces');
        bcnodes = bcnodefacetbl.get('nodes');
        extbcnodefacetbl.faces = repmat(bcfaces, Nd, 1);
        extbcnodefacetbl.nodes = repmat(bcnodes, Nd, 1);
        extbcnodefacetbl = IndexTable(extbcnodefacetbl);
        
        bc.bcnodefacetbl = extbcnodefacetbl;
        bc.linform = vertcat(linforms{:});
        bc.linformvals = vertcat(linformvals{:});
        clear extbcnodefacetbl linforms linformvals

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
        loadstruct.extforce = zeros(tbls.nodefacecoltbl.num, 1);
        clear bc force
        
        assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings);
        clear prop loadstruct
        
        B   = assembly.B  ;
        rhs = assembly.rhs;
        sol = B\rhs;

        % Displacement values at cell centers.
        cellcoltbl = tbls.cellcoltbl;
        n = cellcoltbl.num;

        u = sol(1 : n);
        u = reshape(u, 2, [])';
        
        dnum = u;
        
        % Analytical solution : displacement at cell centers
        dex = NaN(Nc, Nd);
        for idim = 1 : Nd
            dex(:, idim) = u_fun{idim}(cc{:}); 
        end
        
        % Errors in L2 and max - norm
        deL2(iter1) = sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, (dex - dnum).^2)))) / sqrt(sum(sum(bsxfun(@times, G.cells.volumes.^2, ( dex).^2)))); 
        dem(iter1) = max(max(abs(dex - dnum))); 
        
        dostresscomputation = false;
        if dostresscomputation
            san1 = sum([s11(xf( :, 1), xf( :, 2)), s21(xf( :, 1), xf( :, 2))] ...
                       .* G.faces.normals, 2); 
            san2 = sum([s12(xf( :, 1), xf( :, 2)), s22(xf( :, 1), xf( :, 2))] ...
                       .* G.faces.normals, 2)
            stress = md.stress*reshape(dnum', [], 1); 
            
            s_ex = reshape([san1, san2]', [], 1); 
            
            sem(iter1) = max(abs(s_ex - stress))/ max(abs(stress)); 
            fa = reshape(repmat(G.faces.areas, 1, Nd)', [], 1); 
            seL2(iter1) = sqrt(sum(fa.^2 .* (stress - s_ex).^2)) / sqrt(sum(fa.^2 .* s_ex.^2)); 
        end
        
        % Compute VEM solution
        if doVem
            error('bc not implemented yet for VEM');
            [E, nu] = elasticModuloTransform(lambda, mu, 'lam_mu', 'E_nu'); 
            Ev = repmat(E, G.cells.num, 1); 
            nuv = repmat(nu, G.cells.num, 1); 
            C = Enu2C(Ev, nuv, G); 
            % set all boundary to no displacement
            faces = find(any(G.faces.neighbors == 0, 2)); 
            inodes = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 1); 
            nodes = unique(G.faces.nodes(inodes)); 
            el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', zeros(numel(nodes), G.griddim), 'faces', faces,...
                                             'uu_face', zeros(numel(nodes), G.griddim), 'mask', true(numel(nodes), G.griddim)),...
                           'force_bc', []); 
            load = @(coord) - [force_fun{1}(coord( :, 1), coord( :, 2)), force_fun{2}(coord( :, 1), coord( :, 2))]; 
            
            [uVEM, extra] = VEM_linElast(G, C, el_bc, load); 
            
            % Prepare input for analytical functions
            for idim = 1 : Nd
                ncc{idim} = G.nodes.coords(G.cells.nodes, idim);
            end
            % Analytical solution : displacement at ncc
            dnex = NaN(numel(ncc{1}), Nd);
            for idim = 1 : Nd
                dnex(:, idim) = u_fun{idim}(ncc{:}); 
            end        
            w = G.weights.cell_nodes;
            deVEM(iter1) = sqrt(sum(sum(bsxfun(@times, (uVEM(G.cells.nodes, : ) - ...
                                                        dnex).^2, w), 2)))./ ...
                sqrt(sum(sum(bsxfun(@times, dnex.^2, w), 2))); 
            % make CC solution with global interface
            [uu, out] = CC_linElast(G, C, el_bc, load); 
        end
                    
    end

    %% Print convergence rates

    fprintf('Convergence rate for MPSA\n');
    log2(deL2(1 : end - 1)./deL2(2 : end))
    if doVem
        fprintf('Convergence rate for VEM\n');
        log2(deVEM(1 : end - 1)./deVEM(2 : end))
    end
    

    %% Error at last iteration

    fprintf('relative L2 error exact vs mpsa: %g\n', deL2(end));
    if doVem
        fprintf('relative L2 error exact vs vem: %g\n', deVEM(end));
    end

    %% 
    if doVem
        n = 3; 
    else 
        n = 2;
    end
    
    figure(1)
    clf
    set(gcf, 'numbertitle', 'off', 'name', 'DISPLACEMENT')

    subplot(n, 2, 1)
    title('x-mpsa')
    plotCellData(G, dnum(: , 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 2)
    title('y-mpsa')
    plotCellData(G, dnum(: , 2), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 3)
    title('x-exact')
    plotCellData(G, dex(: , 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 4)
    title('y-exact')
    plotCellData(G, dex(: , 2), 'edgecolor', 'none'), colorbar

    if doVem
        subplot(n, 2, 5)
        title('x-vem')
        plotNodeData(G, uVEM(: , 1), 'edgecolor', 'none'), colorbar

        subplot(n, 2, 6)
        title('y-vem')
        plotNodeData(G, uVEM(: , 2), 'edgecolor', 'none'), colorbar
    end
    

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

