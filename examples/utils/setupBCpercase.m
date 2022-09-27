function loadstruct = setupBCpercase(runcase, G, tbls, mappings)
% Boundary conditions

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


    % One linear form per Dirichlet condition

    nodefacetbl        = tbls.nodefacetbl;
    nodefacecoltbl     = tbls.nodefacecoltbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    coltbl             = tbls.coltbl;
    cellcoltbl         = tbls.cellcoltbl;
    
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    d_num = coltbl.num;
    cnfc_num = cellnodefacecoltbl.num;
    
    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
    numnodes = double(diff(G.faces.nodePos));
    numnodes = numnodes(fno);
    facetNormals = G.faces.normals(fno, :);
    facetNormals = bsxfun(@ldivide, numnodes, facetNormals);

    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    facetNormals = sgn.*facetNormals; % Outward normals with respect to cell
                                      % in cellnodefacetbl.
    facetNormals = reshape(facetNormals', [], 1);

    % no volumetric force
    force = zeros(cellcoltbl.num, 1);
    
    
    switch runcase
        
      case {'2d-refinement', '2d-linear', '2d-compaction'}
        extfaces{1} = find(G.faces.centroids(:, 2) == 0);
        n = numel(extfaces{1});
        linforms{1} = repmat([0, 1], n, 1);
        linformvals{1} = zeros(n, 1);

        switch runcase
          case {'2d-linear', '2d-refinement'}
            extfaces{2} = find(G.faces.centroids(:, 1) == 0);
            n = numel(extfaces{2});
            linforms{2} = repmat([1, 0], n, 1);
            linformvals{2} = zeros(n, 1);

          case '2d-compaction'
            extfaces{2} = extfaces{1};
            linforms{2} = repmat([1, 0], n, 1);
            linformvals{2} = zeros(n, 1);

          otherwise
            error('runcase not recognized');
        end

        % Setup force at top, in opposite normal direction
        y = G.faces.centroids(:, 2);
        ymax = max(y);
        extfacetbl.faces = find(y == ymax);
        extfacetbl = IndexArray(extfacetbl);

        [extnodefacetbl, indstruct] = crossIndexArray(nodefacetbl, extfacetbl, {'faces'});
        nodeface_from_extnodeface = indstruct{1}.inds;

        extnodefacecoltbl = crossIndexArray(extnodefacetbl, coltbl, {});

        map = TensorMap();
        map.fromTbl  = cellnodefacecoltbl;
        map.toTbl    = extnodefacecoltbl;
        map.mergefds = {'faces', 'nodes', 'coldim'};

        % Here, we setup a direct construction of extcellnodefacetbl, which allows us
        % to avoid map = map.setup()    
        %    
        %    map = map.setup();
        %    extFacetNormals = map.eval(facetNormals);

        u = (1 : extnodefacetbl.num)';
        v = zeros(nodefacetbl.num, 1);
        v(nodeface_from_extnodeface) = u;
        w = v(nodeface_from_cellnodeface);
        ind = find(w);
        fds = {'cells', 'nodes', 'faces'};
        [~, fdind] = ismember(fds, cellnodefacecoltbl.fdnames);
        extcellnodefacemat = cellnodefacecoltbl.inds(:, fdind);
        extcellnodefacemat = extcellnodefacemat(ind, 1 : 3);

        extcellnodefacetbl = IndexArray([]);
        extcellnodefacetbl = extcellnodefacetbl.setup(fds, extcellnodefacemat);

        extnodeface_from_extcellnodeface = w(ind);
        cellnodeface_from_extcellnodeface = ind;

        extcellnodefacecoltbl = crossIndexArray(extcellnodefacetbl, coltbl, {});

        map.pivottbl = extcellnodefacecoltbl;
        ecnf_num = extcellnodefacetbl.num;
        ecnfc_num = extcellnodefacecoltbl.num;
        [c, i] = ind2sub([d_num, ecnf_num], (1 : ecnfc_num)');
        ind1 = cellnodeface_from_extcellnodeface(i);
        map.dispind1 = sub2ind([d_num, cnfc_num], c, ind1);
        ind2 = extnodeface_from_extcellnodeface(i);
        map.dispind2 = sub2ind([d_num, cnfc_num], c, ind2);
        map.issetup = true;

        extFacetNormals = map.eval(facetNormals);

        map = TensorMap();
        map.fromTbl = extnodefacecoltbl;
        map.toTbl = nodefacecoltbl;
        map.mergefds = {'faces', 'nodes', 'coldim'};

        map = map.setup();

        extforce = map.eval(-extFacetNormals);

        dosourceterm = false;
        if dosourceterm
            error('not done');
            % We setup a source-term
            switch dim
              case 2
                indcell = floor(nx/2) + nx*floor((ny - 1)/2);
                force = [0; 1]; % force in upward direction
              case 3
                indcell = floor(nx/2 + ny/2*nx + nz/2*nx*ny);
                force = [0; 0; 1]; % force in upward direction    
            end

            sourcetbl.cells = indcell;
            sourcetbl.num = numel(indcell);

            sourcetbl = crossIndexArray(sourcetbl, coltbl, {});

            force = tblmap(force, coltbl, sourcetbl, {'coldim'});
            force = tblmap(force, sourcetbl, cellcoltbl, {'cells', 'coldim'});
        end

      case {'3d-linear', '3d-compaction'}
        switch runcase
          case '3d-linear'
            for i = 1 : 3
                extfaces{i} = find(G.faces.centroids(:, i) == 0);
                linform = zeros(3, 1);
                linform(i) = 1;
                n = numel(extfaces{i});
                linforms{i} = repmat(linform, n, 1);
            end
          case '3d-compaction'
            extface = find(G.faces.centroids(:, 3) == 0);
            for i = 1 : 3
                extfaces{i} = extface;
                linform = zeros(3, 1);
                linform(i) = 1;
                n = numel(extfaces{i});
                linforms{i} = repmat(linform, n, 1);
                if i == 1
                    % We impose a translation in the x-direction on the top face
                    linformvals{i} = ones(n, 1);
                else
                    linformvals{i} = ones(n, 1);
                end
            end
          case '3d-compaction-dirichlet'
            extface = find(G.faces.centroids(:, 3) == 0);
            for i = 1 : 3
                extfaces{i} = extface;
                linform = zeros(3, 1);
                linform(i) = 1;
                n = numel(extfaces{i});
                linforms{i} = repmat(linform, n, 1);
                linformvals{i} = zeros(n, 1); 
            end
        end

        % Setup force at top, in opposite normal direction
        y = G.faces.centroids(:, 3);
        ymax = max(y);
        extfacetbl.faces = find(y == ymax);
        extfacetbl = IndexArray(extfacetbl);

        [extnodefacetbl, indstruct] = crossIndexArray(nodefacetbl, extfacetbl, {'faces'});
        nodeface_from_extnodeface = indstruct{1}.inds;

        extnodefacecoltbl = crossIndexArray(extnodefacetbl, coltbl, {});

        map = TensorMap();
        map.fromTbl  = cellnodefacecoltbl;
        map.toTbl    = extnodefacecoltbl;
        map.mergefds = {'faces', 'nodes', 'coldim'};

        map = map.setup();
        extFacetNormals = map.eval(facetNormals);

        map = TensorMap();
        map.fromTbl = extnodefacecoltbl;
        map.toTbl = nodefacecoltbl;
        map.mergefds = {'faces', 'nodes', 'coldim'};

        map = map.setup();
        extforce = map.eval(-extFacetNormals);

      otherwise
        error('runcase not recognized');
    end

    bc.extfaces    = vertcat(extfaces{:});
    bc.linform     = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});

    bc = setupFaceBC(bc, G, tbls);

    loadstruct.bc = bc;
    loadstruct.extforce = extforce;
    loadstruct.force = force;
end
