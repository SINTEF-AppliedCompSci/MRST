function G = createBisectedTetrahedralGrid(Nx,alternate)
    
    [X,Y,Z] = meshgrid(linspace(0,1,Nx(1)+1), linspace(0,1,Nx(2)+1),linspace(0,1,Nx(3)+1));
    p = [X(:), Y(:), Z(:)];
    
    G = mcomputeGeometry(tetrahedralGrid(p));
    
end

