function [G, G_org] = complex3DGrid(opt, grid_case)
%
%
% SYNOPSIS:
%   function [G, G_org] = complex3DGrid(opt, grid_case)
%
% DESCRIPTION:
%
% PARAMETERS:
%
%   opt       -  Parameter structure with fields
%        L            - Physical dimension ('box' and 'grdecl' case)
%        cartDims     - Cartesian dimension ('box' and 'grdecl' case)
%        disturb      - disturbance parameter for the grid ('box' case)
%        triangulate  - If true, the horizontal faces are triangulated.
%        vertical     - Only relevant for Norne grid (straightens up
%                       the pillars, see paper [Andersen et al: http://arxiv.org/abs/1606.09508v1])
%        gtol         - Grid tolerance parameter (used when calling
%                       processGRDECL, see documentation there)
%        ref          - Refinement parameter, only used for Norne grid
%
%   grid_case - 3D Grid cases that are set up here.
%          'box'    - grid set up using squareGrid function
%          'grdecl' - grid set up using Eclipse format
%          'sbed'   - run dataset_bedmodel2() to get a full description
%          'Norne'  - run dataset_norne() to get a full description
%
% RETURNS:
%   G     - Grid that has been created
%   G_org - Original grid without refinement (only 'Norne' case)
%
% EXAMPLE:
%    G=complex3DGrid(struct('vertical',true,'triangulate',false,'ref',1,'gtol',1e-3),'norne');clf,plotGrid(G),view(3)
%    G=complex3DGrid(struct('triangulate',false,'gtol',1e-3),'sbed');clf,plotGrid(G),view(3)
%    G=complex3DGrid(struct('triangulate',true,'gtol',1e-3),'sbed');clf,plotGrid(G),view(3)
%    G=complex3DGrid(struct('triangulate',false,'gtol',1e-3),'sbed');clf,plotGrid(flipGrid(G)),view(3)
%    G=complex3DGrid(struct('cartDims',[3,3,3],'L',[3 3 3],'disturb',0.04,'triangulate',false),'box');clf,plotGrid(G),view(3)
%    G=complex3DGrid([],'grdecl'); clf, plotGrid(G), view(3)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    G_org = [];
    switch grid_case

      case 'box'
        G = squareGrid(opt.cartDims, opt.L, 'grid_type', 'cartgrid', 'disturb', opt.disturb);
        G_old = computeGeometry(G);
        if (opt.triangulate)
            face = unique(G.cells.faces(any(bsxfun(@eq, G.cells.faces(:, 2), [5, 6]), 2), 1));
            G = triangulateFaces(G, [face]');
        end
        G = computeGeometry(G);
        G = createAugmentedGrid(G);

      case 'grdecl'
        griddim=3;
        grdecl = simpleGrdecl([2, 1, 2]*ceil((1e3).^(1/griddim)), 0.15);
        G = processGRDECL(grdecl);
        G = createAugmentedGrid(G);
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

        grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
        grdecl = readGRDECL(grdecl);
        usys   = getUnitSystem('METRIC');
        grdecl = convertInputUnits(grdecl, usys);
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
