function [Gt, rock] = makeTopSurfaceGrid(dim, lengths, depth, poro, perm)
% ----------------------------------------------------------------------------
    Gt   = topSurfaceGrid(computeGeometry(setTopDepth(cartGrid(dim, lengths), depth)));
    rock = averageRock(struct('perm', expand_var(perm, prod(dim)), ...
                              'poro', expand_var(poro, prod(dim))), Gt);
end
