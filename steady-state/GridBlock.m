classdef GridBlock
    %Base class for upscaling a single block

properties
    verbose
    deck      % deck for the block only. May be empty.
    G         % Grid structure for the block. May be periodic.
    rock      % Rock structure for the block. Used for upscaling.
    rockorg   % Original rock structure.
    fluid     % Fluid structure for the block. May be empty.
    periodic  % Boolean. True if block is periodic and false otherwise.
    periodicDims % Dimensions which are made periodic. Defaulted to all.
    bcp       % Periodic boundary conditions. Empty if periodic is false.
    faces     % Side faces
    areas     % Side areas
    lengths   % Vector of characteristic length for each dimension.
    pv        % Pore volume
end

methods
    
    function block = GridBlock(G, rock, varargin)
        block.verbose  = mrstVerbose();
        block.periodic = false;
        block.deck     = [];
        block.faces    = [];
        block.lengths  = [];
        block.areas    = [];
        block.periodicDims = [];
        block = merge_options(block, varargin{:});
        
        block.G       = G; % May be periodic grid
        block.rock    = rock;
        block.rockorg = rock;
        block.pv      = poreVolume(G, rock);
        
        % Set faces if not given
        if isempty(block.faces)
            block.faces = block.findSideFaces(G);
        end
        
        % Compute all face areas
        if isempty(block.areas)
            A = nan(3,2);
            for d=1:3
                for s=1:2
                    A(d,s) = sum(G.faces.areas(block.faces{d}{s}));
                end
            end
            block.areas = A;
        end
        
        % Compute lengths
        if isempty(block.lengths)
            L = nan(3,1);
            for d=1:3
                L(d) = abs(...
                    mean(G.faces.centroids(block.faces{d}{2},d)) - ...
                    mean(G.faces.centroids(block.faces{d}{1},d)) );
            end
            block.lengths = L;
        end
        
        if block.periodic
            % If periodic option is true, then a periodic grid is created
            [block.G, block.bcp] = makePeriodicCartesianGrid(block.G);
        end
        
    end
    
    
    
end


methods (Static)
    
    function faces = findSideFaces(G)
        % Given a grid, find the side faces on each side in each direction.
        % NOTE: This method assumed a high regularity of the grid, and may
        % fail for more complex grids.
        
        bf   = prod(G.faces.neighbors, 2) == 0; % boundary faces
        
        % For each face, find the direction of the largest component of the
        % normal vector. That is, we assume the normal vector maily points
        % in of of the main directions.
        [~, ninx] = max(abs(G.faces.normals), [], 2);
        
        
        faces = cell(3,1);
        ndims = 3;
        for d = 1:3
            odims = [1:d-1 d+1:3]; % other dimensions
            
            fd = find(bf & ninx == d); % boundary face and normal in dir d
            fL = G.faces.centroids(fd, d) < ...
                    mean(G.faces.centroids(fd, d));
            f{1} = fd(fL);  % minimum side ("left" side)
            f{2} = fd(~fL); % maximum side ("right" side)
            
            % We now sort the order of the boundary faces on both sides of
            % the grid, in order for them to "line up" with the
            % corresponding face on the other side of the grid. This, of
            % course, assumes a highly regular grid.
            for i = 1:2 % min, max side ("left" and "right")
                for j = 1:ndims-1 % the other dimensions
                    [~, inx] = sort(G.faces.centroids(f{i}, odims(j)));
                    f{i} = f{i}(inx);
                end
            end
            
            faces{d}    = cell(1,2);
            faces{d}{1} = f{1};
            faces{d}{2} = f{2};
        end
        
    end
    
end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
