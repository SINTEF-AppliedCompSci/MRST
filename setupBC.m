function [D, force] = setupBC(runcase, G, tbls, mappings, varargin)
%% Boundary conditions

    opt = struct('facetNormals', []);
    opt = merge_options(opt, varargin{:}); 
    
    % One linear form per Dirichlet condition

    nodefacetbl        = tbls.nodefacetbl;
    nodefacecoltbl     = tbls.nodefacecoltbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    coltbl             = tbls.coltbl;

    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    d_num = coltbl.num;
    cnfc_num = cellnodefacecoltbl.num;
    
    
    switch runcase
        
      case {'2d-refinement', '2d-linear', '2d-compaction'}
        
        bottomextfaces = find(G.faces.centroids(:, 2) == 0);
        d = [0; 1];
        Dmat{1} = assignLinearForm(bottomextfaces, d, tbls);

        switch runcase
            
          case {'2d-linear', '2d-refinement'}
            leftextfaces = find(G.faces.centroids(:, 1) == 0);
            d = [1; 0];
            Dmat{2} = assignLinearForm(leftextfaces, [1; 0], tbls);
            
          case '2d-compaction'
            d = [1; 0];
            Dmat{2} = assignLinearForm(bottomextfaces, [1; 0], tbls);            

          otherwise
            error('runcase not recognized');
        end

        D = [Dmat{1}, Dmat{2}];

        % Setup force at top, in opposite normal direction
        y = G.faces.centroids(:, 2);
        ymax = max(y);
        extfacetbl.faces = find(y == ymax);
        extfacetbl.num = numel(extfacetbl.faces);

        [extnodefacetbl, indstruct] = crossTable(nodefacetbl, extfacetbl, {'faces'});
        nodeface_from_extnodeface = indstruct{1}.inds;

        extnodefacecoltbl = crossTable(extnodefacetbl, coltbl, {});

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
        cellnodefacemat = convertTableToArray(cellnodefacetbl, fds);
        extcellnodefacemat = cellnodefacemat(ind, 1 : 3);
        extcellnodefacetbl = convertArrayToTable(extcellnodefacemat, fds);
        extnodeface_from_extcellnodeface = w(ind);
        cellnodeface_from_extcellnodeface = ind;

        extcellnodefacecoltbl = crossTable(extcellnodefacetbl, coltbl, {});

        map.pivottbl = extcellnodefacecoltbl;
        ecnf_num = extcellnodefacetbl.num;
        ecnfc_num = extcellnodefacecoltbl.num;
        [c, i] = ind2sub([d_num, ecnf_num], (1 : ecnfc_num)');
        ind1 = cellnodeface_from_extcellnodeface(i);
        map.dispind1 = sub2ind([d_num, cnfc_num], c, ind1);
        ind2 = extnodeface_from_extcellnodeface(i);
        map.dispind2 = sub2ind([d_num, cnfc_num], c, ind2);
        map.issetup = true;

        facetNormals = opt.facetNormals;
        extFacetNormals = map.eval(facetNormals);

        map = TensorMap();
        map.fromTbl = extnodefacecoltbl;
        map.toTbl = nodefacecoltbl;
        map.mergefds = {'faces', 'nodes', 'coldim'};

        map = map.setup();

        force = map.eval(-extFacetNormals);

        dosourceterm = false;
        if dosourceterm
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

            sourcetbl = crossTable(sourcetbl, coltbl, {});

            force = tblmap(force, coltbl, sourcetbl, {'coldim'});
            force = tblmap(force, sourcetbl, cellcoltbl, {'cells', 'coldim'});
        end
      
      case '3d-linear'
        Dmat = cell(3, 1);
        for i = 1 : 3
            extfaces = find(G.faces.centroids(:, i) == 0);
            d = zeros(3, 1);
            d(i) = 1;
            Dmat{i} = assignLinearForm(extfaces, d, tbls);
        end
        D = horzcat(Dmat{:});
        
        % Setup force at top, in opposite normal direction
        y = G.faces.centroids(:, 2);
        ymax = max(y);
        extfacetbl.faces = find(y == ymax);
        extfacetbl.num = numel(extfacetbl.faces);

        [extnodefacetbl, indstruct] = crossTable(nodefacetbl, extfacetbl, {'faces'});
        nodeface_from_extnodeface = indstruct{1}.inds;

        extnodefacecoltbl = crossTable(extnodefacetbl, coltbl, {});

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
        force = map.eval(-extFacetNormals);
      
      otherwise
        error('runcase not recognized');
        
    end
    
end

function D = assignLinearForm(extfaces, linform, tbls)
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;
  
    extfacetbl.faces = extfaces;
    extfacetbl.num = numel(extfacetbl.faces);

    extnodefacetbl = crossTable(nodefacetbl, extfacetbl, {'faces'});

    prod = TensorProd();
    prod.tbl1 = coltbl;
    prod.tbl2 = extnodefacetbl;
    prod.tbl3 = nodefacecoltbl;
    prod = prod.setup();

    D_T = SparseTensor('matlabsparse', true);
    D_T = D_T.setFromTensorProd(linform, prod);
    D = D_T.getMatrix();
    
end
