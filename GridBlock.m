classdef GridBlock
%Base class for upscaling a single block

properties
    verbose
    G         % Grid structure for the block. May be periodic.
    rock      % Rock structure for the block.
    fluid     % Fluid structure for the block. May be empty.
    periodic  % Boolean. True if block is periodic and false otherwise.
    bcp       % Periodic boundary conditions. Empty if periodic is false.
    facesl    % Cell array of left side faces for each dimension.
    facesr    % Cell array of right side faces for each dimension.
    lengths   % Vector of characteristic length for each dimension.
end

methods
    
    function upscaler = GridBlock(G, rock, varargin)
        upscaler.verbose  = mrstVerbose();
        upscaler.periodic = false;
        upscaler.facesl   = [];
        upscaler.facesr   = [];
        upscaler = merge_options(upscaler, varargin{:});
        
        if upscaler.periodic
            % If periodic option is true, then a periodic grid is created
            periodicDims = 1:3; % what edges to make periodic
            [G, upscaler.bcp] = makePeriodicBlock(G, 'dims', periodicDims);
            
            
            
        end
        
        upscaler.G    = G; % May be periodic grid
        upscaler.rock = rock;
        
        % Compute transmissiblity for block
        if upscaler.periodic
            upscaler.T = computeTransGp(G.parent, G, rock);
        else
            upscaler.T = computeTrans(G, rock);
        end
        
        % Compute areas of sides
        areas = nan(ndims,1);
        for i = 1:ndims
            if upscaler.periodic
                areas(i) = sum( G.faces.areas( ...
                    bcp.face(bcp.tags==dims(i))) );
            else
                areas(i) = sum( G.faces.areas(faces{dims(i),2}) );
            end
        end
        
    end
    
    function [data, report] = upscale(upscaler)
        
        report = [];
        wantReport = nargout>1;
        if wantReport
            startTime = tic;
        end
        
        % Perform actual upscaling here
        
        if wantReport
            report.time = toc(startTime);
        end
        
    end
    
    
    function [Gp, bcp] = makePeriodicBlock(G, varargin)
        % Pressure drop in each direction is set to zero
        opt = struct( ...
           'dims',   1:3  ... % dimensions to make periodic
           );
        opt = merge_options(opt, varargin{:});
        
        ndims = numel(opt.dims);
        bcl   = cell(1, ndims); % "left" faces
        bcr   = cell(1, ndims); % "right" faces
        
        for id = 1:ndims
            d = opt.dims(id);
            bcl{id}.face = upscaler.facesl{d};
            bcr{id}.face = upscaler.facesr{d};
        end
        
        % Set pressure drop to zero
        dp = cell(1, ndims);
        dp(:) = {0};
        
        % Create periodic grid
        [Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dp);
    end
    
    function findSideFaces(upscaler)
        
        grid = upscaler.G;
        bf   = prod(grid.faces.neighbors, 2) == 0; % boundary faces
        
        % For each face, find the direction of the largest component of the
        % normal vector. That is, we assume the normal vector maily points
        % in of of the main directions.
        [~, ninx] = max(abs(upscaler.G.faces.normals), [], 2);
        
        
        fleft = cell(1,3);
        fright = cell(1,3);
        for d = 1:3
            odims = [1:d-1 d+1:3]; % other dimensions
            
            fd = find(bf & ninx == d); % boundary face and normal in dir d
            fL = grid.faces.centroids(fd, d) < ...
                    mean(grid.faces.centroids(fd, d));
            f{1} = fd(fL);  % minimum side ("left" side)
            f{2} = fd(~fL); % maximum side ("right" side)
            
            % We now sort the order of the boundary faces on both sides of
            % the grid, in order for them to "line up" with the
            % corresponding face on the other side of the grid. This, of
            % course, assumes a highly regular grid.
            for i = 1:2 % min, max side ("left" and "right")
                for j = 1:ndims-1 % the other dimensions
                    [~, inx] = sort(grid.faces.centroids(f{i}, odims(j)));
                    f{i} = f{i}(inx);
                end
            end
            
            fleft{d} = f{1};
            fright{d} = f{2};
        end
        
        upscaler.facesl = fleft;
        upscaler.facesr = fright;
        
    end
    
    
end
    
end



