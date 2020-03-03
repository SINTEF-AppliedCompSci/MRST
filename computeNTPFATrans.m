function [M, M_noflow] = computeNTPFATrans(G, rock)
    
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
            locneig = G.faces.neighbors(i, 1);
            f = OSflux{i, j};
            if  ~isempty(f)
                
                ncells   = size(f, 1) - 1;
                locfaces = repmat(i, ncells, 1);
                loccells2 = repmat(locneig, ncells, 1);
                loccells1 = f(1 : end - 1, 1); 
                locvals   = f(1 : end - 1, 2);
                
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

    vals = [vals{1}; vals{2}];

    celltbl.cells = (1 : nc)';
    celltbl = IndexTable(celltbl);

    facetbl.faces = (1 : nf)';
    facetbl = IndexTable(facetbl);

    
    prod = TensorProd();
    prod.tbl1 = cell2facetbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = facetbl;
    prod.replacefds2 = {{'cells', 'cells1'}};
    prod.reducefds = {'cells1'};

    prod = prod.setup();

    tens = SparseTensor();
    tens = tens.setFromTensorProd(vals, prod);
    
    M = tens.getMatrix();

    isBoundary = any(G.faces.neighbors == 0, 2);
    intfacetbl.faces = find(~isBoundary);
    intfacetbl = IndexTable(intfacetbl);

    map = TensorMap();
    map.fromTbl = facetbl;
    map.toTbl = intfacetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    tens2 = SparseTensor();
    tens2 = tens2.setFromTensorMap(map);
    
    tens = tens2*tens;

    M_noflow = tens.getMatrix();
    
end
