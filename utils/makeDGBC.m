function bc = makeDGBC(disc, bc)

    G = disc.G;
    assert(all(strcmpi(bc.type, 'flux')));
    [faces, cells] = boundaryFaces(G);
    
    ix = any(faces == bc.face',2);
    cells = [ones(nnz(ix),1),cells(ix)];
    
    checkUpstr = false(numel(bc.face), 1);
    
    bc.cell = cells;
    bc.checkUpstr = checkUpstr;
    st = struct('sdof', bc.sat, ...
                'degree', zeros(numel(bc.face),1));
    dd = disc;
    dd.G.cells.num = numel(bc.face);
    st = dd.updateDisc(st);
            
    bc.state = st;

end