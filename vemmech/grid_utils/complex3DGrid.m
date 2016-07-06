function [G, G_org] = complex3DGrid(opt, grid_case)

    G_org = [];
    switch grid_case

      case 'box'
        G = squareGrid(opt.cartDims, opt.L, 'grid_type', 'cartgrid', 'disturb', opt.disturb);
        G_old = computeGeometry(G);
        if(opt.triangulate)
            face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [5, 6]), 2), 1));
            G = triangulateFaces(G, [face]');
        end
        G = computeGeometry(G);
        G = mrstGridWithFullMappings(G);

      case 'grdecl'

        grdecl = simpleGrdecl([2, 1, 2]*ceil((1e3).^(1/griddim)), 0.15);
        G = processGRDECL(grdecl);
        G = mrstGridWithFullMappings(G);
        G = computeGeometry(G);

      case 'sbed'
        % Load the grid for the sbed model as eclipse input
        grdecl = readGRDECL(fullfile(getDatasetPath('BedModel2'), ...
                                     'BedModel2.grdecl'));
        grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));

        grdecl.ZCORN = min(grdecl.ZCORN, 5)-1;
        grdecl.ZCORN = max(grdecl.ZCORN, 1)-1;
        grdecl_c = cutGrdecl(grdecl, [1 15;1 15; 1 333]); % consider only a
                                                          % part of the sbed model
        G = processGRDECL(grdecl_c, 'Tolerance', opt.gtol);
        if (opt.triangulate)
            faces = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [5, 6]), 2), 1));
            G = triangulateFaces(G, faces');
        end
        figure();
        plotGrid(G);

      case 'norne'
        % Load the grid for the Norne model as eclipse input
        grdecl = readGRDECL(fullfile(getDatasetPath('BedModel2'), ...
                                     'BedModel2.grdecl'));
        grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
        grdecl = cutGrdecl(grdecl, [10 25;35 55;1 22]);
        if (opt.vertical)
            grdecl_org = verticalGrdecl(grdecl);
        else
            grdecl_org = grdecl;
        end
        G_org = processGRDECL(grdecl_org);
        if (opt.triangulate)
            faces = unique(G_org.cells.faces(any(bsxfun(@eq, G_org.cells.faces(:, 2), [5, 6]), 2), 1));
            G_org = triangulateFaces(G_org, faces');
        end
        grdecl = padGrdecl(grdecl, [true, true, true], [60 50; 40 40; 10 10]*3, 'relative', true);
        if (opt.vertical)
            grdecl = verticalGrdecl(grdecl);
        end
        grdecl = refineGrdeclLayers(grdecl, [1 1], opt.ref);
        grdecl = refineGrdeclLayers(grdecl, [grdecl.cartDims(3) grdecl.cartDims(3)], opt.ref);
        G = processGRDECL(grdecl, 'Tolerance', opt.gtol);
        G = computeGeometry(G);
        if(opt.triangulate)
            faces = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [5, 6]), 2), 1));
            G = triangulateFaces(G, faces');
        end
        figure();
        clf;
        plotGrid(G);

      otherwise
        error();
    end
end