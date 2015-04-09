classdef Upscaler
%Base class for upscaling classes

properties
    verbose
    G
    rock
    fluid
    periodic
end

methods
    
    function upscaler = Upscaler(G, rock, fluid, varargin)
        upscaler.verbose  = mrstVerbose();
        upscaler.periodic = false;
        upscaler = merge_options(upscaler, varargin{:});
        
        upscaler.G     = G;
        upscaler.rock  = rock;
        upscaler.fluid = fluid;
    end
    
    function [data, report] = upscale(upscaler, partition)
        
        nblocks = max(partition);
        assert(numel(partition)==upscaler.G.cells.num, ...
            'Invalid partition');
        
        report = [];
        wantReport = nargout>1;
        if wantReport
            report.time = nan(nblocks,1);
        end
        
        %data = struct([]);
        
        % Loop over blocks and perform upscaling on each block
        for b=1:nblocks
            
            if upscaler.verbose
                fprintf('Block number %d of %d\n', b, nblocks);
            end
            
            if wantReport
                t = tic;
            end
            
            % Create grid, rock and fluid for sub block
            cells = find(partition==b);
            [BG, BRock, BFluid] = createBlock(upscaler, cells);
            
            % Perform upscaling
            data(b) = upscaler.upscaleBlock(BG, BRock, BFluid); %#ok<AGROW>
            
            if wantReport
                report.time(b) = toc(t);
            end
        end
        
    end
    
    function [BG, BRock, BFluid] = createBlock(upscaler, cells)
        % Create grid, rock and fluid for the sub block represented by the
        % given cells.
        
        % Create block grid
        BG = extractSubgrid(upscaler.G, cells);
        %BG.cartDims = subGriddim(cells);
        %BG.cells.indexMap = (1 : BG.cells.num)';
        
        % Create periodic grid
        if upscaler.periodic
            error('Periodic grid not implemented')
        end
        
        % Create block rock
        BRock.perm = upscaler.rock.perm(cells, :);
        BRock.poro = upscaler.rock.poro(cells, :);
        if isfield(upscaler.rock, 'cr')
            BRock.cr = upscaler.rock.cr;
        end
        
        % Create block fluid
        BFluid = createBlockFluid(upscaler.fluid, cells);
    end
    
    function data = upscaleBlock(upscaler, BG, BRock, BFluid)
        error('Method needs to be overridden');
    end
    
end
    
end

