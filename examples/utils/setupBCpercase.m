function loadstruct = setupBCpercase(runcase, G, tbls, mappings, extras, varargin)
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

    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    % One linear form per Dirichlet condition

    nodefacetbl        = tbls.nodefacetbl;
    nodefacevectbl     = tbls.nodefacevectbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacevectbl = tbls.cellnodefacevectbl;
    vectbl             = tbls.vectbl;
    cellvectbl         = tbls.cellvectbl;
    celltbl            = tbls.celltbl;
    
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    d_num = vectbl.num;
    cnfc_num = cellnodefacevectbl.num;
    
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

    % force and extforce will be set to zero if nothing is set otherwise
    force    = [];
    extforce = [];
    
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

        extnodefacevectbl = crossIndexArray(extnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

        [extcellnodefacetbl, indstruct] = crossIndexArray(extnodefacetbl, cellnodefacetbl, {'nodes', 'faces'});

        extnodeface_from_extcellnodeface  = indstruct{1}.inds;
        cellnodeface_from_extcellnodeface = indstruct{2}.inds;

        extcellnodefacevectbl = crossIndexArray(extcellnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

        map = TensorMap();
        map.fromTbl  = cellnodefacevectbl;
        map.toTbl    = extnodefacevectbl;
        map.mergefds = {'faces', 'nodes', 'vec'};
        
        map.pivottbl = extcellnodefacevectbl;
        ecnf_num = extcellnodefacetbl.num;
        ecnfc_num = extcellnodefacevectbl.num;
        [c, i] = ind2sub([d_num, ecnf_num], (1 : ecnfc_num)');
        ind1 = cellnodeface_from_extcellnodeface(i);
        map.dispind1 = sub2ind([d_num, cnfc_num], c, ind1);
        ind2 = extnodeface_from_extcellnodeface(i);
        map.dispind2 = sub2ind([d_num, cnfc_num], c, ind2);
        map.issetup = true;

        extFacetNormals = map.eval(facetNormals);

        map = TensorMap();
        map.fromTbl  = extnodefacevectbl;
        map.toTbl    = nodefacevectbl;
        map.mergefds = {'faces', 'nodes', 'vec'};

        if useVirtual
            
            map.pivottbl = extnodefacevectbl;

            [vec, i] = ind2sub([d_num, extnodefacetbl.num], (1 : extnodefacevectbl.num)');
            map.dispind1 = (1 : extnodefacevectbl.num)';
            map.dispind2 = sub2ind([d_num, nodefacetbl.num], vec, nodeface_from_extnodeface(i));
            map.issetup = true;
            
        else
            
            map = map.setup();
            
        end

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
            sourcetbl       = IndexArray(sourcetbl);

            map = TensorMap();
            map.fromTbl  = celltbl;
            map.toTbl    = sourcetbl;
            map.mergefds = {'cells'};

            cell_from_source = map.getDispatchInd();
            
            sourcevectbl = crossIndexArray(sourcetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

            map = TensorMap();
            map.fromTbl  = vectbl;
            map.toTbl    = sourcevectbl;
            map.mergefds = {'vec'};

            if useVirtual
                
                map.pivottbl = sourcevectbl;
                
                [vec, i] = ind2sub([vectbl.num, sourcetbl.num], (1 : sourcevectbl.num)');
                map.dispind1 = vec;
                map.dispind2 = (1 : sourcevectbl.num)';

                map.issetup = true;
            else
                
                map = map.setup();
                
            end
            
            force = map.eval(force);

            map = TensorMap();
            map.fromTbl  = sourcevectbl;
            map.toTbl    = cellvectbl;
            map.mergefds = {'cells', 'vec'};
            
            if useVirtual
                
                map.pivottbl = sourcevectbl;
                
                [vec, i] = ind2sub([vectbl.num, sourcetbl.num], (1 : sourcevectbl.num)');
                map.dispind1 = (1 : sourcevectbl.num);
                map.dispind2 = sub2ind([vectbl.num, celltbl.num], vec, cell_from_source(i));;

                map.issetup = true;
            else
                
                map = map.setup();
                
            end
            
            force = tblmap(force, sourcetbl, cellvectbl, {'cells', 'vec'});
            
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

        extnodefacevectbl = crossIndexArray(extnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

        map = TensorMap();
        map.fromTbl  = cellnodefacevectbl;
        map.toTbl    = extnodefacevectbl;
        map.mergefds = {'faces', 'nodes', 'vec'};

        if useVirtual

            [extcellnodefacetbl, indstruct] = crossIndexArray(extnodefacetbl, cellnodefacetbl, {'nodes', 'faces'});
            
            extnodeface_from_extcellnodeface  = indstruct{1}.inds;
            cellnodeface_from_extcellnodeface = indstruct{2}.inds;
            
            extcellnodefacevectbl = crossIndexArray(extcellnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

            map.pivottbl = extcellnodefacevectbl;

            [vec, i] = ind2sub([vectbl.num, extcellnodefacetbl.num], (1 : extcellnodefacevectbl.num)');
            map.dispind1 = sub2ind([vectbl.num, cellnodefacetbl.num], vec, cellnodeface_from_extcellnodeface(i));
            map.dispind2 = sub2ind([vectbl.num, extnodefacetbl.num], vec, extnodeface_from_extcellnodeface(i));

            map.issetup = true;
            
        else
            
            map = map.setup();
            
        end
        extFacetNormals = map.eval(facetNormals);

        map = TensorMap();
        map.fromTbl  = extnodefacevectbl;
        map.toTbl    = nodefacevectbl;
        map.mergefds = {'faces', 'nodes', 'vec'};

        if useVirtual
            
            map.pivottbl = extnodefacevectbl;

            [vec, i] = ind2sub([vectbl.num, extnodefacetbl.num], (1 : extnodefacevectbl.num)');
            map.dispind1 = (1 : extnodefacevectbl.num)';
            map.dispind2 = sub2ind([vectbl.num, nodefacetbl.num], vec, nodeface_from_extnodeface(i));

            map.issetup = true;
        else

            map = map.setup();
            
        end
        
        extforce = map.eval(-extFacetNormals);

      case '3d-gravity'

        ind = 1;
        for i = 1 : 3
            for j = 1 : 2
                addface = true;
                if j == 1
                    if i == 3
                        % The bottom face is free. (Note that for geological reservoir, bottom face corresponds to top face.)
                        addface = false;
                    else
                        v = min(G.faces.centroids(:, i));
                    end
                else
                    v = max(G.faces.centroids(:, i));
                end

                if addface
                    
                    extfaces{ind} = find(abs(G.faces.centroids(:, i) - v) < 1e-5);
                    
                    n = numel(extfaces{ind});

                    linform = zeros(3, 1);
                    linform(i) = 1;
                    linforms{ind}    = repmat(linform, n, 1);

                    linformvals{ind} = zeros(n, 1);

                    ind = ind + 1;
                    
                end
            end
        end

        % for ind = 1 : numel(extfaces)
        %     figure
        %     plotGrid(G, 'facevecor', 'none');
        %     plotFaces(G, extfaces{ind});
        % end
        
        vols = G.cells.volumes;

        ztbl.vec = 3;
        ztbl = IndexArray(ztbl);

        [cellztbl, indstruct] = crossIndexArray(celltbl, ztbl, {});
        cell_from_cellztbl = indstruct{1}.inds;
        vec_from_cellztbl  = indstruct{2}.inds;
        
        map = TensorMap();
        map.fromTbl  = cellztbl;
        map.toTbl    = cellvectbl;
        map.mergefds = {'cells', 'vec'};

        if useVirtual
            
            map.pivottbl = cellztbl;

            map.dispind1 = (1 : cellztbl.num)'
            map.dispind2 = sub2ind([vectbl.num, celltbl.num], ...
                                   vec_from_cellztbl, ...
                                   cell_from_cellztbl);
            
            map.issetup = true;
            
        else
            
            map = map.setup();
            
        end
        rho = extras.rho;
        force = map.eval(rho*vols*10);
        
      otherwise
        
        error('runcase not recognized');
        
    end

    bc.extfaces    = vertcat(extfaces{:});
    bc.linform     = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});

    bc = setupFaceBC(bc, G, tbls);

    if isempty(force)
        force = zeros(cellvectbl.num, 1);
    end

    if isempty(extforce)
        extforce = zeros(nodefacevectbl.num, 1);
    end
    
    loadstruct.bc       = bc;
    loadstruct.extforce = extforce;
    loadstruct.force    = force;
    
end
