function G = squareGrid(cartDims, L, varargin)
%
% SYNOPSIS:
%   G = squareGrid(cartDims, L, varargin)
%
% DESCRIPTON:
% Make square with different irregular grid types starting from cartesian.
%
% run script exploreSquareGrid.m to see the different possibilities
%
% REQUIRED PARAMETERS:
%   cartDims  - Cartesian dimensions
%   L         - Physical lengths
%
% OPTIONAL PARAMETERS:
%
%   grid_type - Possible types are  'cartgrid' 'triangle' 'pebi' 'boxes[1-4]'
%               'mixed[1-4]'. Run the scripts given in the comments below to
%               see how the grid look likes. Only the 'cartgrid' type is
%               implemented for 3D. The others apply uniquely for 2D.
%   disturb   - The grids can be twisted, essentially to break the Cartesian
%               structure. The parameter disturb determines the amount of
%               twisting that is introduced.
%
% RETURNS:
%
%   G   - Grid structure

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


    opt = struct('grid_type' , 'cartgrid', ...
                 'disturb'   , 0.0);
    opt = merge_options(opt, varargin{:});
    opt.L = L;
    opt.cartDims = cartDims;

    if(numel(opt.cartDims) == 2)
        [nx, ny] = deal(opt.cartDims(1), opt.cartDims(2));
        [Lx, Ly] = deal(opt.L(1), opt.L(2));
    else
        assert(strcmp(opt.grid_type, 'cartgrid'))
    end

    switch opt.grid_type

      case 'cartgrid'
        % Standard Cartesian grid
        G = cartGrid(opt.cartDims, opt.L);
        if(opt.disturb > 0)
            G = twister(G, opt.disturb);
        end

      case 'triangle'
        % Grid made of triangles
        G_tmp = cartGrid(opt.cartDims, opt.L);
        if(opt.disturb > 0)
            G_tmp = twister(G_tmp, opt.disturb);
        end
        p = [G_tmp.nodes.coords(:, 1), G_tmp.nodes.coords(:, 2)];
        t = delaunayn(p);
        G = triangleGrid(p, t);

      case 'pebi'
        % Pebi grid, dual of delaunay grid.
        G_tmp = cartGrid(opt.cartDims, opt.L);
        if(opt.disturb > 0)
            G_tmp = twister(G_tmp, opt.disturb);
        end
        p = [G_tmp.nodes.coords(:, 1), G_tmp.nodes.coords(:, 2)];
        t = delaunayn(p);
        G = triangleGrid(p, t);
        G = pebi(G);
        G = removeShortEdges(G, min(opt.L)/(max(opt.cartDims)*10));
        G = sortEdges(G);

      case   {'boxes1', 'boxes2', 'boxes3', 'boxes4'}
        % boxes side by side
        % 4 alternatives
        % Run :
        %
        % figure()
        % plotGrid(squareGrid([4, 4],[3, 3], 'grid_type', 'boxed4', 'disturb', 0.1))
        %
        % and change 'boxes4' to to the other parameters to explore the
        % different alternatives that are proposed.
        G1 = cartGrid(opt.cartDims, opt.L);
        G2 = cartGrid(2*opt.cartDims+1, opt.L);
        G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
        if(str2num(opt.grid_type(end))>1)
            G3 = cartGrid(opt.cartDims+2, opt.L);
            G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
            if(str2num(opt.grid_type(end))>2)
                G4 = cartGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)]);
                dx = 0.05;
                G5 = cartGrid([opt.cartDims(1)*5, 1], [3*opt.L(1), dx]);
                G  = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
                G  = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
                G  = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
                if(str2num(opt.grid_type(end))>3)
                    G6 = cartGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx]);
                    G = glue2DGrid(G, translateGrid(G6, [- opt.L(1), -opt.L(2)]));
                    G = glue2DGrid(G, translateGrid(G6, [3*opt.L(1), -opt.L(2)]));
                end
            end
        end
        G.cells = rmfield(G.cells, 'indexMap');
        if(opt.disturb > 0)
            G = twister(G, opt.disturb);
        end
        G.nodes.coords(:, 1) = G.nodes.coords(:, 1)-min(G.nodes.coords(:, 1));
        G.nodes.coords(:, 2) = G.nodes.coords(:, 2)-min(G.nodes.coords(:, 2));
        Lx = max(G.nodes.coords(:, 1));
        Ly = max(G.nodes.coords(:, 2));
        G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*opt.L(1)./Lx;
        G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
        G = sortEdges(G);

      case   {'mixed1', 'mixed2', 'mixed3', 'mixed4'}
        % Mixes cells with different types, Cartesian, pebi and triangles
        % 4 alternatives
        % Run :
        %
        % figure()
        % plotGrid(squareGrid([4, 4],[3, 3], 'grid_type', 'mixed4', 'disturb', 0.1))
        %
        % and change 'mixed4' to to the other parameters to explore the
        % different alternatives that are proposed.
        G1 = squareGrid(opt.cartDims, opt.L, 'grid_type', 'cartgrid', 'disturb', opt.disturb);
        G2 = squareGrid(opt.cartDims+1, opt.L, 'grid_type', 'pebi', 'disturb', opt.disturb);
        G = glue2DGrid(G1, translateGrid(G2, [opt.L(1), -0.0]));
        if(str2num(opt.grid_type(end))>1)
            G3 = squareGrid(opt.cartDims+2, opt.L, 'grid_type', 'triangle', 'disturb', 0.02);
            G = glue2DGrid(G, translateGrid(G3, 2*[opt.L(1), 0.0]));
            if(str2num(opt.grid_type(end))>2)
                G4 = cartGrid(opt.cartDims+[7, 1], [3*opt.L(1), opt.L(2)]);
                dx = 0.05;
                G5 = cartGrid([opt.cartDims(1)*5, 1], [3*opt.L(1), dx]);
                G = glue2DGrid(G, translateGrid(G5, [0, opt.L(2)]));
                G = glue2DGrid(G, translateGrid(G4, [0, opt.L(2)+dx]));
                G = glue2DGrid(G, translateGrid(G4, [0, -opt.L(2)]));
                if(str2num(opt.grid_type(end))>3)
                    G7 = squareGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx], 'grid_type', 'triangle', 'disturb', opt.disturb);
                    G6 = squareGrid([3, 10], [opt.L(1), 3*opt.L(2)+dx], 'grid_type', 'pebi', 'disturb', opt.disturb);
                    G = glue2DGrid(G, translateGrid(G6, [- opt.L(1), -opt.L(2)]));
                    G = glue2DGrid(G, translateGrid(G7, [3*opt.L(1), -opt.L(2)]));
                end
            end
        end
        G.cells = rmfield(G.cells, 'indexMap');
        if(opt.disturb > 0)
            G = twister(G, opt.disturb);
        end
        G.nodes.coords(:, 1) = G.nodes.coords(:, 1)-min(G.nodes.coords(:, 1));
        G.nodes.coords(:, 2) = G.nodes.coords(:, 2)-min(G.nodes.coords(:, 2));
        Lx = max(G.nodes.coords(:, 1));
        Ly = max(G.nodes.coords(:, 2));
        G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*opt.L(1)./Lx;
        G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*opt.L(2)./Ly;
        G = sortEdges(G);

      otherwise
        error('no such grid type')
    end

    if(~strcmp(opt.grid_type, 'cartgrid'))
        Lx  = max(G.nodes.coords(:, 1));
        Ly  = max(G.nodes.coords(:, 2));
        tag = zeros(G.faces.num, 1);
        tol = 1e-3/max(opt.cartDims)*min(opt.L);
        G   = computeGeometry(G);
        tag(abs(G.faces.centroids(:, 1))<tol)    = 1;
        tag(abs(G.faces.centroids(:, 1)-Lx)<tol) = 2;
        tag(abs(G.faces.centroids(:, 2))<tol)    = 3;
        tag(abs(G.faces.centroids(:, 2)-Ly)<tol) = 4;
        if(G.griddim == 3)
            Lz = max(G.nodes.coords(:, 3));
            tag(abs(G.faces.centroids(:, 3))<tol)    = 3;
            tag(abs(G.faces.centroids(:, 3)-Ly)<tol) = 4;
        end
        G.cells.faces = [G.cells.faces(:, 1), tag(G.cells.faces(:, 1))];
    end
    G = createAugmentedGrid(G);

end