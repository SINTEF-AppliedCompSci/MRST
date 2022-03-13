classdef Upscaler
    %Base class for upscaling classes

properties
    
    verbose
    timeest % Estimate time remaining
    
    deck
    G
    rock
    fluid
    periodic
    periodicDims
    partition
    blocks
    blockMap % Index map for the blocks. Equal index means equal block.
    
end

methods
    
    function upscaler = Upscaler(G, rock, varargin)
        upscaler.verbose   = mrstVerbose();
        upscaler.periodic  = false;
        upscaler.periodicDims = [];
        upscaler.deck      = [];
        upscaler.fluid     = [];
        upscaler.partition = [];
        upscaler.blocks    = [];
        upscaler.blockMap  = [];
        upscaler = merge_options(upscaler, varargin{:});
        
        upscaler.G     = G;
        upscaler.rock  = rock;
    end
    
    function [blockdata, globdata, report] = upscale(upscaler)
        % Given some partition, we loop over the coarse blocks in the grid
        % and upscale each of them in turn.
        
        globdata = [];
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
                'coarse blocks.\n'], nBlocksReq, nBlocksTot);
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
                fprintf(['Will only upscale %d of %d coarse blocks as '...
                    'a block map is given.\n'], nBlocksUp, nBlocksReq);
            elseif nBlocksReq>1
                fprintf(['No block map is given, so all blocks '...
                    'will be upscaled.\n']);
            end
        end
        
        % Store report
        if wantReport
            report.numBlocksTot = nBlocksTot;
            report.reqBlockInx = requestedBlocks;
            report.upsBlockInx = blocksToUpscale;
            report.blocks = cell(nBlocksUp,1);
        end
        
        if upscaler.verbose
            % Line break before starting
            fprintf('\n');
        end
        
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
            [blockdata(i), blockReport] = upscaler.upscaleBlock(block); %#ok<AGROW>
            blockTime = toc(t);
            
            % Store report
            if wantReport
                blockReport.block.celldim = block.G.cartDims;
                blockReport.block.physdim = block.lengths(:)';
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
        
        % Check for values outside range
        if wantReport && upscaler.verbose
            try % Use try block to be sure we don't break something here
                outsideRange = cellfun(@(x) x.relperm.valsOutsideRange, ...
                    report.blocks);
                outsideRange = sum(outsideRange);
                if outsideRange > 0
                    fprintf(['There are %d blocks where the upscaled '...
                        'relperm had values outside valid range!\n'], ...
                        outsideRange);
                end
            catch
            end
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
        
        
        if ~isempty(upscaler.deck)
            
            % We have a deck available, and so we create the block using
            % the deck directly instead of using the grid and fluid.
            
            blockDeck = Upscaler.createBlockDeck(upscaler.deck, cells);
            
            b = initEclipseGrid(blockDeck);
            try
                b = mcomputeGeometry(b); % Use MEX version if possible
            catch %#ok<CTCH>
                b = computeGeometry(b); % Fallback; use Matlab version
            end
            r = initEclipseRock(blockDeck);
            
            f = [];
            if ~isempty(upscaler.fluid)
                f = initDeckADIFluid(blockDeck);
            end
            
            % Create block object
            block = GridBlock(b, r, 'fluid', f, 'deck', blockDeck, ...
                'periodic', upscaler.periodic, ...
                'periodicDims', upscaler.periodicDims);
            
        else
            
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
                f = createBlockFluid(upscaler.fluid, cells);
            end
            
            % Create block object
            block = GridBlock(b, r, 'fluid', f, ...
                'periodic', upscaler.periodic, ...
                'periodicDims', upscaler.periodicDims);
            
        end
        
    end
    
    function [data, report] = upscaleBlock(upscaler, block)
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
    
    function db = createBlockDeck(deck, cells)
        % Create a sub deck for the block based on the given cells
        % NOTE: This function assumes a Cartesian grid.
        
        if isempty(deck)
            db = [];
            return;
        end
        
        if isfield(deck.RUNSPEC, 'cartDims')
            cartDims = deck.RUNSPEC.cartDims;
        elseif isfield(deck.RUNSPEC, 'DIMENS')
            cartDims = deck.RUNSPEC.DIMENS;
        else
            error('Unable to find dimensions of grid from deck');    
        end
        
        db = deck;
        nc = prod(cartDims); % number of cells in whole grid
        
        % Compute block dimensions
        [I,J,K]   = ind2sub(cartDims, cells);
        dim       = @(V) max(V)-min(V)+1;
        bGridDims = [dim(I) dim(J) dim(K)];
        
        % RUNSPEC
        if isfield(deck.RUNSPEC, 'cartDims')
            db.RUNSPEC.cartDims = bGridDims;
        end
        if isfield(deck.RUNSPEC, 'DIMENS')
            db.RUNSPEC.DIMENS = bGridDims;
        end
        
        % GRID
        assert(all(isfield(deck.GRID, {'DX','DY','DZ'})), ...
            'This function only accepts grids given by DX, DY, DZ.');
        cfns = {'DX','DY','DZ','PERMX','PERMY','PERMZ', 'PORO'};
        for i=1:numel(cfns)
            fn = cfns{i};
            if isfield(deck.GRID, fn) && numel(deck.GRID.(fn))==nc
                db.GRID.(fn) = deck.GRID.(fn)(cells);
            end
        end
        if isfield(db.GRID, 'TOPS')
            % We specifically remove TOPS. This is not needed in the
            % upscaling, but it may cause errors if transported directly to
            % a subgrid.
            db.GRID = rmfield(db.GRID, 'TOPS');
        end
        if isfield(deck.GRID, 'cartDims')
            db.GRID.cartDims = bGridDims;
        end
        
        % PROPS
        % No nothing
        
        % REGIONS
        cfns = {'PVTNUM','SATNUM'};
        for i=1:numel(cfns)
            fn = cfns{i};
            if isfield(deck.REGIONS, fn) && numel(deck.REGIONS.(fn))==nc
                db.REGIONS.(fn) = deck.REGIONS.(fn)(cells);
            end
        end
        
        % SOLUTION
        % No nothing
        
        % SUMMARY
        % No nothing
        
        % SCHEDULEf
        % No nothing
        
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
