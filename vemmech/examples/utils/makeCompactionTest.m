function [el_bc, load] = makeCompactionTest(G, opt, varargin)
%
%
% SYNOPSIS:
%   function [el_bc, load] = makeCompactionTest(G, opt)
%
% DESCRIPTION: Set up the compaction test by computing the elastic boundary
% conditions and the load term.
%
% PARAMETERS:
%   G   - Grid
%   opt - Structure with fields:
%          'gravity_load' : Includes gravity in load
%          'hanging'      : no displacement on the sides
%          'free_top'     : no force applied on top
%          'top_load'     : constant pressure applied on top
%          'islinear'     : impose some linear displacement
%                           (in vertical direction)
%
% OPTIONAL PARAMETERS:
%  'gravity'   - value for gravity (default gravity = 10)
%  'density'   - value for density in kg/m^3 (default density = 3000)
%  'top_force' - value for top force (default top_force = 30000)
%          
% RETURNS:
%   el_bc - Elastic boundary condition structure. It contains the fields
%             'disp_bc'  : displacement boundary condition. It contains the 
%                          fields
%                  'nodes'    : nodes where the displacement condition is applied   
%                  'uu'       : value for the displacement
%
%                  The two following fields are not used in the VEM implementation
%                  but becomes relevant for other methods such as MPSA, see paper
%
%                  'faces'    : faces displacement
%                  'uu_faces' : value for the displacement
%
%                  'mask'     : if false then displacement values that are
%                               imposed in given Cartesian directions are in
%                               fact ignored.
%             'force_bc'  : force boundary condition applied on faces. It contains the 
%                           fields
%                  'faces' : faces where the force is applied
%                  'force' : value of the force that is applied
%    
%   load  - load function, can be evaluated anywhere, independently of the
%           grid structure.

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

    opt2 = struct('gravity'  , 10, ... 
                  'density'  , 3000, ...
                  'top_force', 30000,...
                  'rolling_vertical', false,...
                  'tol',sqrt(eps));
    opt2 = merge_options(opt2, varargin{:});

    
    sides = {'XMin', 'XMax', 'YMin', 'YMax', 'ZMin', 'ZMax'};
    for j = 1 : G.griddim
        Lmax = max(G.faces.centroids(:, j));
        Lmin = min(G.faces.centroids(:, j));
        x = [Lmin, Lmax];
        for i = 1 : 2
             if(~isfield(G, 'cartDims'))
                faces = find(abs(G.faces.centroids(:, j)-x(i))<opt2.tol);
                assert(all(any(G.faces.neighbors(faces,:)==0,2)));
             else
                mside         = sides{2*(j-1)+i};
                tmp           = pside([], G, mside, 0);
                faces         = tmp.face;
             end
             bc{i + (j - 1)*2} = addBC([], faces, 'pressure', 0);
        end
    end

    % Find the node of the different sides and set up the elastisity boundary conditions

    for i = 1 : 2*G.griddim
        inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face+1)-1);
        nodes = unique(G.faces.nodes(inodes));
        disp_bc = struct('nodes'  , nodes, ...
                         'uu'     , 0, ...
                         'faces'  , bc{i}.face, ...
                         'uu_face', 0, ...
                         'mask'   , true(numel(nodes), G.griddim));
        bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
    end

    % Set up the gravity load
    
    density = opt2.density; 
    grav = opt2.gravity;    
    gdir = zeros(1, G.griddim);
    gdir(end) = 1;
    if (opt.islinear)
        origo  = mean(G.cells.centroids, 1);
        fac    = (max(G.cells.centroids(:, G.griddim))-min(G.cells.centroids(:, G.griddim)))*1000;
        bcdisp = @(x) bsxfun(@minus, bsxfun(@times, x, [0 0 1]), origo)./fac;
    else
        bcdisp = @(x)(x*0.0);  
    end
    if (opt.gravity_load)
        load = @(x)(-(grav*density)*repmat(gdir, size(x, 1), 1));
    else    
        load = @(x)(-(0*density)*repmat(gdir, size(x, 1), 1));
    end

    % Set the Dirichlet boundary conditions at the selected sides

    bc_el_sides = bc;
    if (~opt.hanging || opt.islinear)
        for i = 1 : 2*(G.griddim - 1)
            % On vertical boundary faces : Rolling in the tangential directions, that is we
            % only impose zero displacement in the normal direction.
            bc_el_sides{i}.el_bc.disp_bc.mask(:, G.griddim) = false;
            if (G.griddim == 3 && ~opt2.rolling_vertical)
                if (i <= 2)
                    bc_el_sides{i}.el_bc.disp_bc.mask(:, 2) = false;
                elseif (i <= 4)
                    bc_el_sides{i}.el_bc.disp_bc.mask(:, 1) = false;
                end
            end
        end
        for i = 2*G.griddim
            % On the bottom boundary faces : Rolling in horizontal direction,
            % no displacement in vertical direction.
            bc_el_sides{i}.el_bc.disp_bc.mask(:, 1 : (G.griddim - 1)) = false;
        end
        for i = 2*G.griddim - 1
            % On the top boundary faces : 
            if (opt.free_top)
                bc_el_sides{i} = []; 
            else    
                bc_el_sides{i}.el_bc.disp_bc.mask(:, 1 : (G.griddim - 1)) = false; 
            end   
        end
    else
        for i = 2*G.griddim-2 : 2*G.griddim
            bc_el_sides{i} = [];
        end  
    end
    % Collect the boundary conditions
    nodes = [];
    faces = [];
    mask  = [];

    for i = 1 : numel(bc)
        if (~isempty(bc_el_sides{i}))
            nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes];
            faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces];
            mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask];  
        end
    end

    disp_node = bcdisp(G.nodes.coords(nodes, :));
    disp_faces = bcdisp(G.faces.centroids(faces, :));
    disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask);

    if (opt.top_load)
        fvec = zeros(1, G.griddim);
        fvec(G.griddim) = 1;
        H = max(G.nodes.coords(:, ...
                               G.griddim))-min(G.nodes.coords(:, G.griddim));
        top_force = opt2.top_force;
        face_force = @(x) H*top_force*repmat(-fvec, size(x, 1), 1);
        faces = bc{2*G.griddim-1}.face;
        % Make force boundary structure. NB: force unit is  Pa/m^3.
        force_bc = struct('faces', faces, 'force', face_force(G.faces.centroids(faces, :)));
    else
        force_bc = []; 
    end

    % Final structure fo boundary conditions
    el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);

end
