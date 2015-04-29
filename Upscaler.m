classdef Upscaler
%Base class for upscaling classes

properties
    
    verbose
    timeest % Estimate time remaining
    
    G
    rock
    fluid
    periodic
    partition
    blocks
    blockMap % Index map for the blocks. Equal index means equal block.
    
    cellinx % Boolean. True if ALL fluid properties takes 'cellInx'.
end

methods
    
    function upscaler = Upscaler(G, rock, varargin)
        upscaler.verbose   = mrstVerbose();
        upscaler.periodic  = false;
        upscaler.fluid     = [];
        upscaler.partition = [];
        upscaler.blocks    = [];
        upscaler.blockMap  = [];
        upscaler.cellinx   = false;
        upscaler = merge_options(upscaler, varargin{:});
        
        upscaler.G     = G;
        upscaler.rock  = rock;
    end
    
    function [data, report] = upscale(upscaler)
        % Given some partition, we loop over the coarse blocks in the grid
        % and upscale each of them in turn.
        
        startTime = tic;
        wantReport = nargout > 1;
        
        % The partition gives the coarse block number for each fine cell
        p  = upscaler.partition;
        nBlocksTot = numel(unique(p)); % Total number of coarse blocks
        
        % Check given partition
        assert(~isempty(p), 'Empty partition');
        assert(numel(p)==upscaler.G.cells.num, 'Invalid partition');
        assert(max(p)==nBlocksTot, 'Invalid partition numbering');
        
        if upscaler.verbose
            fprintf(['Partition divides %d cells into %d coarse '...
                'blocks.\n'], upscaler.G.cells.num, nBlocksTot);
        end
        
        % Block numbers requested to be upscaled (if user do not wish to
        % upscale all coarse blocks in partition).
        requestedBlocks = upscaler.blocks;
        if isempty(requestedBlocks)
            requestedBlocks = (1:nBlocksTot)'; % all blocks
        end
        nBlocksReq = numel(requestedBlocks); % Total number of blocks
        
        if upscaler.verbose && nBlocksReq<nBlocksTot
            fprintf(['Only requested to upscale %d of totally %d '...
                'coarse blocks.\n'], upscaler.G.cells.num, nBlocksTot);
        end
        
        % If given, the blockMap gives identical block, which may reduce
        % the number of blocks needed to be upscaled
        if ~isempty(upscaler.blockMap)
            bm = upscaler.blockMap(requestedBlocks);
            [~, inx] = unique(bm); % Choose one cell from each block
            blocksToUpscale = requestedBlocks(inx);
        else
            blocksToUpscale = requestedBlocks;
        end
        nBlocksUp = numel(blocksToUpscale);
        
        if upscaler.verbose
            if nBlocksUp<nBlocksReq
                fprintf(['Need only upscale %d of %d coarse blocks as '...
                    'some blocks are idential.\n'], nBlocksUp, nBlocksReq);
            else
                fprintf(['No blocks are identical, so all blocks '...
                    'need to be upscaled.\n']);
            end
        end
        
        % Store report
        if wantReport
            report.numBlocksTot = nBlocksTot;
            report.reqBlockInx = requestedBlocks;
            report.upsBlockInx = blocksToUpscale;
            report.blocks = cell(nBlocksUp,1);
        end
        
        % Line break before starting
        fprintf('\n');
        
        % Loop over blocks and perform upscaling on each block
        for i = 1:nBlocksUp
            
            t = tic;
            
            if upscaler.verbose
                fprintf('Block number %d of %d\n', i, nBlocksUp);
            end
            
            b = blocksToUpscale(i); % Current block
            
            % Create grid, rock and fluid for sub block
            cells = find(p==b);
            block = upscaler.createBlock(cells); %#ok<FNDSB>
            
            setupTime = toc(t);
            if upscaler.verbose
                fprintf('  Setup block:  %6.3fs\n', setupTime);
            end
            
            % Perform upscaling
            [data(i), blockReport] = upscaler.upscaleBlock(block); %#ok<AGROW>
            blockTime = toc(t);
            
            % Store report
            if wantReport
                report.blocks{i} = blockReport;
            end
            
            % Print info
            if upscaler.verbose
                fprintf('  Total time:   %6.3fs\n', blockTime);
                if upscaler.timeest
                    totalTime  = toc(startTime);
                    estTimeRem = totalTime*(nBlocksUp/i - 1);
                    timeStr = Upscaler.timingString(estTimeRem);
                    fprintf('  Estimated time left is %s\n', timeStr);
                end
                fprintf('\n');
            end
        end
        
        totalTime = toc(startTime);
        
        % Store report
        if wantReport
            report.time = totalTime;
        end
        
        if upscaler.verbose
            timeStr = Upscaler.timingString(totalTime);
            fprintf('Completed upscaling of %d blocks in %s.\n', ...
                nBlocksUp, timeStr);
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


methods (Static)
    
    function str = timingString(seconds)
        h = floor(seconds/(60*60)); seconds = seconds - h*60*60;
        m = floor(seconds/60);      seconds = seconds - m*60;
        s = floor(seconds);
        if h>0
            str = sprintf('%dh %dm %ds', h, m, s);
        elseif m>0
            str = sprintf('%dm %ds', m, s);
        else
            str = sprintf('%ds',s);
        end
    end
    
end


    
end

