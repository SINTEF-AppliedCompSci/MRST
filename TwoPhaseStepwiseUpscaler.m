classdef TwoPhaseStepwiseUpscaler < OnePhaseUpscaler
%Two phase upscaling

properties
    method1
    method2
    dim1
    dim2
    nrelperm
    pcow
    npcow
    
    saveStep1
end

methods
    
    function upscaler = TwoPhaseStepwiseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid);
        
        upscaler.method1   = [];
        upscaler.method2   = [];
        upscaler.dim1      = [];
        upscaler.dim2      = [];
        upscaler.nrelperm  = 20;
        upscaler.pcow      = true; % Upscale pcow or not
        upscaler.npcow     = 100;
        
        upscaler.saveStep1 = false; % save data from first upscaling
        
        upscaler = merge_options(upscaler, varargin{:});
        
        % The absolute permeability will be evaluated in the second
        % dimension chosen.
        % upscaler.dims = upscaler.dim2;
        
    end
       
    function [data, report] = upscaleBlock(upscaler, block)
        
        wantReport = nargout > 1;
        startTime  = tic;
        
        % Separate special upscaling in the z-direction is implemented, but
        % is set as unreachable code for now. The method is questionable.
        upscalez   = true;
        
        assert(~isempty(upscaler.method1) && ...
            ~isempty(upscaler.method2), 'Methods must be given.');
        assert(~isempty(upscaler.dim1) && ...
            ~isempty(upscaler.dim2), 'Dimensions must be given.');
        
        % Perform one phase upscaling first
        [data, rep] = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        if wantReport
            report.onephase = rep;
        end
        
        %------------------------------------------------------------------
        % Relative permeability upscaling, part 1
        % (Upscale in direction 1, using method 1)
        %------------------------------------------------------------------
        
        % Direction 1 is typically z
        
        t = tic;
        t_relperm = t;
        ijk = gridLogicalIndices(block.G);
        p1  = ijk{upscaler.dim2};
        
        if upscalez
            dims1 = [upscaler.dim2, upscaler.dim1];
        else
            dims1 = upscaler.dim2;
        end
        
        upscaler1 = TwoPhaseUpscaler(block.G, block.rock, block.fluid, ...
            'partition', p1, 'method', upscaler.method1, ...
            'dims', dims1, 'nrelperm', upscaler.nrelperm, ...
            'pcow', true, 'verbose', false);
        if wantReport
            [updata1, rep] = upscaler1.upscale();
            report.relperm.method = 'stepwise';
            report.relperm.step1  = rep;
        else
            updata1 = upscaler1.upscale();
        end
        relperm1Time = toc(t);
        
        if upscaler.verbose
            fprintf('  Rel.perm #1:  %6.3fs\n', relperm1Time);
        end
        
        %------------------------------------------------------------------
        % Relative permeability upscaling, part 2
        % (Upscale in direction 2, using method 2)
        %------------------------------------------------------------------
        
        t = tic;
        
        % Create coarse grid from partition1
        CG = generateCoarseGrid(block.G, p1);
        CG = coarsenGeometry(CG);
        CG = addNodeDataToCoarseGrid(CG); % needed in extractSubgrid
        cdims = block.G.cartDims;
        cdims(upscaler.dim1) = 1;
        CG.cartDims = cdims;
        CG.cells.indexMap = (1:CG.cells.num)';
        
        % Create coarse rock and fluid
        perm       = reshape([updata1.perm], numel(dims1), [])';
        CRock.perm = perm(:,1); % use perm from dim2 pressure drop
        CRock.poro = reshape([updata1.poro],[],1);
        warning('off', 'SteadyState:MultidimRelperm');
        CFluid = createUpscaledFluid(block.fluid, updata1, p1);
        warning('on', 'SteadyState:MultidimRelperm');
        
        % Upscale in second direction
        p2 = ones(CG.cells.num,1);
        upscaler2 = TwoPhaseUpscaler(CG, CRock, CFluid, ...
            'partition', p2, 'method', upscaler.method2, ...
            'dims', upscaler.dim2, 'nrelperm', upscaler.nrelperm, ...
            'pcow', false, 'verbose', false, 'cellinx', true);
        if wantReport
            [updata2, rep] = upscaler2.upscale();
            report.relperm.step2 = rep;
        else
            updata2 = upscaler2.upscale();
        end
        krO = updata2.krO;
        krW = updata2.krW;
        
        if upscalez
            % Upscale relperm in z-direction in a special manner
            sW   = krO{1}(:,1);
            krzW = nan(numel(sW), 1);
            krzO = nan(numel(sW), 1);
            for is=1:numel(sW)
                % Get relperm value in each intermediate block for current sW
                int  = @(kr) interp1(kr(:,1), kr(:,2), sW(is));
                intW = arrayfun(@(x) int(extendTab(x.krW{2})), updata1);
                intO = arrayfun(@(x) int(extendTab(x.krO{2})), updata1);

                % The upscaled relperm value is taken as the arithmetic mean
                krzW(is) = mean(intW);
                krzO(is) = mean(intO);
            end
            
            % Enforce endpoints to zero (this is not true in general when
            % simply using the arithmetic mean)
            krzW(1)   = 0;
            krzO(end) = 0;
            
            % Using x-direction also in y-direcion (or vica versa)
            krW = {krW{1}, krW{1}, [sW krzW]};
            krO = {krO{1}, krO{1}, [sW krzO]};
        end
        
        % Save to structure
        data.krO = krO;
        data.krW = krW;
        relperm2Time = toc(t);
        
        if upscaler.saveStep1
            % Save data from the first upscaling step, if requested
            data.krStep1 = updata1;
            data.krStep1.method    = upscaler.method1;
            data.krStep1.dim       = upscaler.dim1;
            data.krStep1.partition = p1;
        end
        
        if upscaler.verbose
            fprintf('  Rel.perm #2:  %6.3fs\n', relperm2Time);
            totalTime = toc(t_relperm);
            fprintf('  Rel.perm tot: %6.3fs\n', totalTime);
        end
        
        
        %------------------------------------------------------------------
        % Capillary pressure upscaling
        %------------------------------------------------------------------
        if upscaler.pcow
            [data, rep] = upPcOW(block, data, ...
                'npointsmax', upscaler.npcow);
            if wantReport
                report.pcow = rep;
            end
            if upscaler.verbose
                fprintf('  Cap.pres:     %6.3fs\n', rep.time);
            end
        end
        
        if wantReport
            totalTime   = toc(startTime);
            report.time = totalTime;
        end
        
        
    end
    
end

end


