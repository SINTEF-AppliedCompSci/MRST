classdef Upscaler
%Base class for upscaling classes

properties
    verbose
    G
    rock
    fluid
    periodic
    partition
    blocks
    cellinx % fluid property
end

methods
    
    function upscaler = Upscaler(G, rock, varargin)
        upscaler.verbose   = mrstVerbose();
        upscaler.periodic  = false;
        upscaler.fluid     = [];
        upscaler.partition = [];
        upscaler.blocks    = [];
        upscaler.cellinx   = false;
        upscaler = merge_options(upscaler, varargin{:});
        
        upscaler.G     = G;
        upscaler.rock  = rock;
    end
    
    function data = upscale(upscaler)
        
        startTime = tic;
        
        p  = upscaler.partition;
        bl = upscaler.blocks;
        if isempty(bl)
            bl = unique(p);
        end
        nblocks = numel(bl);
        assert(numel(p)==upscaler.G.cells.num, ...
            'Invalid partition');
        
        % Loop over blocks and perform upscaling on each block
        for i = 1:nblocks
            
            b = bl(i); % Current block
            
            if upscaler.verbose
                fprintf('Block number %d of %d\n', b, nblocks);
            end
            
            % Create grid, rock and fluid for sub block
            t = tic;
            cells = find(p==b);
            block = upscaler.createBlock(cells); %#ok<FNDSB>
            t = toc(t);
            if upscaler.verbose
                fprintf('  Setup block:  %6.3f sec.\n', t);
            end
            
            % Perform upscaling
            t = tic;
            data(i) = upscaler.upscaleBlock(block); %#ok<AGROW>
            if upscaler.verbose
                t = toc(t);
                fprintf('  Total time:   %6.3f sec.\n', t);
            end
        end
        
        if upscaler.verbose
            totalTime = toc(startTime);
            fprintf('Completed upscaling of %d blocks in %1.2f sec.\n', ...
                nblocks, totalTime);
        end
        
    end
    
    function block = createBlock(upscaler, cells)
        % Create grid, rock and fluid for the sub block represented by the
        % given cells.
        
        % Create block grid
        b = extractSubgrid(upscaler.G, cells);
        ijk = gridLogicalIndices(upscaler.G);
        b.cartDims = cellfun(@(x) max(x(cells))-min(x(cells))+1, ijk);
        b.cells.indexMap = (1 : numel(cells))';
        
        % Create block rock
        r.perm = upscaler.rock.perm(cells, :);
        r.poro = upscaler.rock.poro(cells, :);
        if isfield(upscaler.rock, 'cr')
            r.cr = upscaler.rock.cr;
        end
        
        % Create block fluid
        f = [];
        if ~isempty(upscaler.fluid)
            f = createBlockFluid(upscaler.fluid, cells, ...
                'cellInx', upscaler.cellinx);
        end
        
        % Create block object
        block = GridBlock(b, r, 'fluid', f);
    end
    
    function data = upscaleBlock(upscaler, block)
        error('Method needs to be overridden');
    end
    
end


    
end

