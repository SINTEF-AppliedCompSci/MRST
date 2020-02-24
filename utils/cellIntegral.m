function int_u = cellIntegral(u, cubature, cells)
    % Integrate integrand u over cells using a given cubature
    % int_u(i) = (int_{cell(i)} u dv)/|cell(i)|
    if nargin < 3 || isinf(cells)
        % Empty cells means all cells in grid
        cells = (1:cubature.G.cells.num)';
    end
    % Get cubature for all cells, transform coordinates to ref space
    W = cubature.getCubature(cells, 'cell');
    int_u = W*u;
end