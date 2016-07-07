function [div] = VEM_div(G)
    if (G.griddim == 3)
        div = VEM3D_div(G);
    else
        assert(G.griddim == 2)
        div = VEM2D_div(G);
    end
end
