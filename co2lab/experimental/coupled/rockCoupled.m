function rockCoupled = rockCoupled(G_coupled, rock, rock2D, region3D)
    % Take rock and rock2D and combine into a suitable rock for a combined
    % 2d/3d grid.
    rockCoupled.perm = rock.perm(G_coupled.cells.indexMap,:);
    rockCoupled.poro = rock.poro(G_coupled.cells.indexMap,:);
    rockCoupled.perm(~region3D) = rock2D.perm(G_coupled.cells.mapTopSurface(~region3D));
end
