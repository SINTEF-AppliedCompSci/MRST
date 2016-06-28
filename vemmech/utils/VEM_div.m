function [div] = VEM_div(G)
% For now define all the primary operators with out units i.e without
% volum, area, length

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
    if (G.griddim == 3)
        div = VEM3D_div(G);
    else
        assert(G.griddim == 2)
        div = VEM2D_div(G);
    end
end
