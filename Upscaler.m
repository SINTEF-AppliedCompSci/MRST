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
end

methods
    
    function upscaler = Upscaler(G, rock, varargin)
        upscaler.verbose   = mrstVerbose();
        upscaler.periodic  = false;
        upscaler.fluid     = [];
        upscaler.partition = [];
        upscaler.blocks    = [];
        upscaler = merge_options(upscaler, varargin{:});
        
        upscaler.G     = G;
        upscaler.rock  = rock;
    end
    
    function [data, report] = upscale(upscaler)
        
        p  = upscaler.partition;
        bl = upscaler.blocks;
        if isempty(bl)
            bl = unique(p);
        end
        nblocks = numel(bl);
        assert(numel(p)==upscaler.G.cells.num, ...
            'Invalid partition');
        
        report = [];
        wantReport = nargout>1;
        if wantReport
            report.time = nan(nblocks,1);
        end
        
        % Loop over blocks and perform upscaling on each block
        for i = 1:nblocks
            
            b = bl(i); % Current block
            
            if upscaler.verbose
                fprintf('Block number %d of %d\n', b, nblocks);
            end
            
            if wantReport
                t = tic;
            end
            
            % Create grid, rock and fluid for sub block
            cells = find(p==b);
            block = createBlock(upscaler, cells);
            
            % Perform upscaling
            data(i) = upscaler.upscaleBlock(block); %#ok<AGROW>
            
            if wantReport
                report.time(i) = toc(t);
            end
        end
        
    end
    
    function block = createBlock(upscaler, cells)
        % Create grid, rock and fluid for the sub block represented by the
        % given cells.
        
        % Create block grid
        b = extractSubgrid(upscaler.G, cells);
        
        % Create block rock
        r.perm = upscaler.rock.perm(cells, :);
        r.poro = upscaler.rock.poro(cells, :);
        if isfield(upscaler.rock, 'cr')
            r.cr = upscaler.rock.cr;
        end
        
        % Create block fluid
        f = [];
        if ~isempty(upscaler.fluid)
            f = createBlockFluid(upscaler.fluid, cells);
        end
        
        % Create block object
        block = GridBlock(b, r, 'fluid', f);
    end
    
    function data = upscaleBlock(upscaler, block)
        error('Method needs to be overridden');
    end
    
end
    
end

