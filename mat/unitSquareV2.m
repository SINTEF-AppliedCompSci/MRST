function G = unitSquareV2(nx, ny)

G = cartGrid([nx, ny], [1, 1]);

X = G.nodes.coords;
internalNodes = find((sum((X == 0) + (X == 1),2) == 0));
Nni = numel(internalNodes);
G.nodes.coords(internalNodes,:) = G.nodes.coords(internalNodes,:) + ...
                     [random('Normal', 0, 1/(8*(nx-1)), Nni, 1),    ...
                      random('Normal', 0, 1/(8*(ny-1)), Nni, 1)];
                      
end