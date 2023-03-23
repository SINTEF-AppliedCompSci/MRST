function W_ghost = addGhostWellsEverywhere(G, W, rock)
    ijk = gridLogicalIndices(G);
    c = [W.cells];
    ij_w = [ijk{1}(c) ijk{2}(c)];

    W_ghost = [];
    for i = 1:5:G.cartDims(1)
        for j = 1:5:G.cartDims(2)
            if ismember([i, j], ij_w, 'rows')
                continue
            end
            W_ghost = verticalWell(W_ghost, G, rock, i, j, [], 'Val', eps, 'Type', 'rate', 'Name', 'MrGhost');
        end
    end
end
