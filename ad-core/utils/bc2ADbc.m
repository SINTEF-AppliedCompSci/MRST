function bcAD=bc2ADbc(G,bc)
    bcAD=bc;
    assert(all(any(G.faces.neighbors(bc.face,:)==0,2))); %bc schould be on a boundary
    bc_cell=sum(G.faces.neighbors(bc.face,:),2);    
    if(false)
        % this used full matrix dimensions
        bcAD.cell2bcface=sparse(bc.face,bc_cell,1,G.faces.num,G.cells.num);
        bcAD.bcface2cell=bcAD.cell2bcface';
    else
       % this uses matrixed reduced to size of number of faces, but keep
       % the dimension fore cells. one could use smaller dimesion for cells
       %  but then the matrixes would not be transpose but one had to
       %  introduce a matrix for cells to bc_cells.??
       nbc=numel(bc.face);
       bc_face_index=[1:nbc]';
       bcAD.cell2bcface=sparse(bc_face_index,bc_cell,1,nbc,G.cells.num);
       bcAD.bcface2cell=bcAD.cell2bcface';         
    end
end