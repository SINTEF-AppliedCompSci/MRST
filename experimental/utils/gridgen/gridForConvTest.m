function G = gridForConvTest(Nx,gridType)
    %{
      Grid type:
      1: Cartesian
      2: Triangles by alternating bisection of triangles
      3: Equilateral triangles
      4: Triangles by uniform bisection
      5: Tetrehedral grid  
    %}
    
    Nd = numel(Nx);

    switch gridType
      case 1
        G = computeGeometry(cartGrid(Nx,ones(1,Nd)));
      case 2
        G = createBisectedTriangleGrid(Nx,1);
      case 3
        G = createEquilateralTriangleGrid(Nx);
      case 4    
        G = createBisectedTriangleGrid(Nx,0);
      case 5
        assert(Nd == 3)
        G = createBisectedTetrahedralGrid(Nx);
      otherwise
        error('gridType not recognized');
    end
end

