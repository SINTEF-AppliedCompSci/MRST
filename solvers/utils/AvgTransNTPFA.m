function T = AvgTransNTPFA(G, OSflux)
    
    nf = G.faces.num;
    nc = G.cells.num;

    celltbl.cells = (1 : nc)';
    celltbl = IndexArray(celltbl);
    
    facetbl.faces = (1 : nf)';
    facetbl = IndexArray(facetbl);
    
    internal = (1 : nf)';
    internal(~all(G.faces.neighbors ~= 0, 2)) = [];
    intfacetbl.faces = internal;
    intfacetbl = IndexArray(intfacetbl);
    
    T = cell(2, 1);
    
    %% We unpack th OSflux to multi-index format using IndexArray
    cells   = cell(2, 1);
    faces   = cell(2, 1);
    weights = cell(2, 1);
    for j = 1 : 2
        cells{j} = [];
        faces{j} = [];
        weights{j} = [];
    end
    
    for i = 1 : nf
        
        for j = 1 : 2
            
            d = OSflux{i, j};
            d = d(:, 1 : end);
            
            if ~isempty(d)
                faces{j} = [faces{j}; repmat(i, size(d, 1), 1)];
                cells{j} = [cells{j}; d(:, 1)];
                weights{j} = [weights{j}; d(:, 2)];
            end
            
        end
        
    end
    
    lincellfacetbls = cell(2, 1);
    
    for j = 1 : 2
        
        clear lincellfacetbl
        lincellfacetbl.cells = cells{j};
        lincellfacetbl.faces = faces{j};
        lincellfacetbl = IndexArray(lincellfacetbl);
        
        % We only consider the internal faces
        linintcellfacetbl = crossIndexArray(lincellfacetbl, intfacetbl , {'faces'});
        
        map = TensorMap();
        map.fromTbl = lincellfacetbl;
        map.toTbl = linintcellfacetbl;
        map.mergefds = {'cells', 'faces'};
        map = map.setup();
        
        weights{j} = map.eval(weights{j});
        lincellfacetbls{j} = linintcellfacetbl;
        
    end
    
    %% We have to change the weights of the "inward" cell to opposite sign
    
    intfaces = intfacetbl.get('faces');
    
    for j = 1 : 2
        
        clear cellfacetbl;
        cellfacetbl.faces = intfaces;
        cellfacetbl.cells = G.faces.neighbors(intfaces, j);
        cellfacetbl = IndexArray(cellfacetbl);

        map = TensorMap();
        map.fromTbl = cellfacetbl;
        map.toTbl = lincellfacetbls{j};
        map.mergefds = {'cells', 'faces'};
        map = map.setup();
        
        chsign = logical(map.eval(ones(cellfacetbl.num, 1)));
        
        weights{j}(chsign) = -weights{j}(chsign);
        
    end        
    
    %% We setup the flux mappings
    
    for j = 1 : 2
        
        prod = TensorProd();
        prod.tbl1 = lincellfacetbls{j};
        prod.tbl2 = celltbl;
        prod.tbl3 = intfacetbl;
        prod.reducefds = {'cells'};
        prod = prod.setup();
        
        T{j} = SparseTensor();
        T{j} = T{j}.setFromTensorProd(weights{j}, prod);
        
        T{j} = T{j}.getMatrix();
        
    end
    
end

