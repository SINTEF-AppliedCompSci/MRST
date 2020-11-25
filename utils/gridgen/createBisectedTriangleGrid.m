function G = createBisectedTriangleGrid(Nx, alternate)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    Nd = numel(Nx); 
    Nc = prod(Nx); % Number of Cartesian blocks, half the number of actual triangles

    [X, Y] = meshgrid(linspace(0, 1, Nx(1) + 1), linspace(0, 1, Nx(2) + 1)); 
    X = X'; 
    Y = Y'; 
    p = [X(:), Y(:)]; 

    tmp = (1:Nx(1))'; 
    i1 = tmp; 
    i2 = tmp + 1; 
    i3 = Nx(1) + 2 + tmp; 
    i4 = Nx(1) + 1 + tmp; 

    t1 = reshape([i1 i2 i3, i1 i3 i4]', Nd + 1, [])'; 

    if alternate
        swap = zeros(Nc, 1); 
        swap(2:2:end) = 1; 
    else
        swap = zeros(Nc, 1); 
    end

    swap = zeros(Nx(1), 1); 
    if alternate
        swap(1:2:Nx(1)) = 1; 
    end

    t = t1; 
    for iter1 = 1 : Nx(2) - 1
        t = [t; t1 + iter1 * (Nx(1) + 1)]; 
        
        if alternate
            s = zeros(Nx(1), 1); 
            if mod(iter1, 2) == 0
                s(1:2:end) = 1; 
            else
                s(2:2:end) = 1; 
            end
            swap = [swap; s]; 
        end
        
    end

    if any(swap)
        to = t; 
        for iter1 = 1 : Nc
            if swap(iter1)
                t(iter1 * Nd - 1, :) = t(iter1 * Nd - 1, :) - [0 0 1]; 
                t(iter1 * Nd, :) = t(iter1 * Nd, :) + [1 0 0]; 
                a = []; 
            end
        end
    end

    g = triangleGrid(p, t); 
    G = computeGeometry(g); 
end
