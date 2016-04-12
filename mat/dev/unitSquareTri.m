function G = unitSquareTri(cartDims, gridLim)

dx = gridLim(1)/cartDims(1); dy = gridLim(1)/cartDims(2);
[x, y] = meshgrid(0:dx:gridLim(1), 0:dy:gridLim(2));

P = [x(:), y(:)];
G = triangleGrid(P);

end