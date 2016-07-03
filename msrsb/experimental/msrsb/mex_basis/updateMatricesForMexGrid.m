function g = updateMatricesForMexGrid(A, CG, g)
    N = CG.cells.num;
    for i = 1:N
        c = g.cells{i} + 1;
        g.subsys{i} = A(c, c);
    end
end