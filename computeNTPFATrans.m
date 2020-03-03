function [fluxop, fluxop_neumann, A] = computeNTPFATrans(G, rock)
    
    bc.face = [];
    bc_nfvm = convertBC2FlowNTPFA(G, bc);
    
    interpFace = findHAP(G, rock, bc_nfvm);
    OSflux = findOSflux(G, rock, bc_nfvm, interpFace);
    
    nf = G.faces.num;
    nc = G.cells.num;

    cells1 = cell(2, 1);
    cells2 = cell(2, 1);
    faces = cell(2, 1);
    vals  = cell(2, 1);

    for j = 1 : 2
        faces{j} = [];
        cells1{j} = [];
        cells2{j} = [];
        vals{j}  = [];
        for i = 1 : nf
            locneig = G.faces.neighbors(i, j);
            f = OSflux{i, j};
            if  ~isempty(f)
                ncells   = size(f, 1) - 1;
                locfaces = repmat(i, ncells, 1);
                loccells2 = repmat(locneig, ncells, 1);
                loccells1 = f(1 : end - 1, 1); 
                locvals   = f(1 : end - 1, 2);
                
                locvals(loccells2 == loccells1) = - locvals(loccells2 == loccells1);
                
                faces{j}  = vertcat(faces{j}, locfaces);
                cells1{j} = vertcat(cells1{j}, loccells1);
                cells2{j} = vertcat(cells2{j}, loccells2);
                vals{j}   = vertcat(vals{j}, locvals);
                
            end
        end
    end

    clear cell2facetbl
    cell2facetbl.cells1 = [cells1{1}; cells1{2}];
    cell2facetbl.cells2 = [cells2{1}; cells2{2}];
    cell2facetbl.faces = [faces{1}; faces{2}];
    cell2facetbl = IndexTable(cell2facetbl);

    vals = -[0.5*vals{1}; -0.5*vals{2}];

    celltbl.cells = (1 : nc)';
    celltbl = IndexTable(celltbl);

    facetbl.faces = (1 : nf)';
    facetbl = IndexTable(facetbl);
    
    multicellfacetbl = projTable(cell2facetbl, {'cells1', 'faces'});
    
    map = TensorMap();
    map.fromTbl = cell2facetbl;
    map.toTbl = multicellfacetbl;
    map.mergefds = {'cells1', 'faces'};
    map = map.setup();
    
    vals = map.eval(vals);

    prod = TensorProd();
    prod.tbl1 = multicellfacetbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = facetbl;
    prod.replacefds1 = {{'cells1', 'cells'}};
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    fluxT = SparseTensor();
    fluxT = fluxT.setFromTensorProd(vals, prod);
    
    fluxop = fluxT.getMatrix();

    N=G.faces.neighbors;
    intInx = all(N ~= 0, 2);
    intfacetbl.faces = find(intInx);
    intfacetbl = IndexTable(intfacetbl);

    map = TensorMap();
    map.fromTbl = facetbl;
    map.toTbl = intfacetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    intT = SparseTensor();
    intT = intT.setFromTensorMap(map);
    
    fluxT_neumann = intT*fluxT;

    fluxop_neumann = fluxT_neumann.getMatrix();
    
    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl = IndexTable(cellfacetbl);
    
    
    cno = cellfacetbl.get('cells');
    fno = cellfacetbl.get('faces');
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;
    
    prod = TensorProd();
    prod.tbl1 = cellfacetbl;
    prod.tbl2 = facetbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'faces'};
    prod = prod.setup();
    
    divT = SparseTensor();
    divT = divT.setFromTensorProd(sgn, prod);
    
    AT = divT*fluxT;
    
    A = AT.getMatrix();
    
    
end
