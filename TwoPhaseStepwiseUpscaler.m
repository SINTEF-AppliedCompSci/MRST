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
        
    end
       
    function data = upscaleBlock(upscaler, block)
        
        assert(~isempty(upscaler.method1) && ...
            ~isempty(upscaler.method2), 'Methods must be given.');
        assert(~isempty(upscaler.dim1) && ...
            ~isempty(upscaler.dim2), 'Dimensions must be given.');
        
        % Perform one phase upscaling first
        data = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        
        %------------------------------------------------------------------
        % Relative permeability upscaling, part 1
        % (Upscale in direction 1, using method 1)
        %------------------------------------------------------------------
        t = tic;
        t_relperm = t;
        ijk = gridLogicalIndices(block.G);
        p1  = ijk{upscaler.dim2};
        upscaler1 = TwoPhaseUpscaler(block.G, block.rock, block.fluid, ...
            'partition', p1, 'method', upscaler.method1, ...
            'dims', upscaler.dim1, 'nrelperm', upscaler.nrelperm, ...
            'pcow', true, 'verbose', false);
        updata1 = upscaler1.upscale();
        
        if upscaler.verbose
            t = toc(t);
            fprintf('  Rel.perm #1:  %6.3f sec.\n', t);
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
        CRock.perm = reshape([updata1.perm],[],1);
        CRock.poro = reshape([updata1.poro],[],1);
        CFluid = createUpscaledFluid(block.fluid, updata1, p1);
        
        % Upscale in second direction
        p2 = ones(CG.cells.num,1);
        upscaler2 = TwoPhaseUpscaler(CG, CRock, CFluid, ...
            'partition', p2, 'method', upscaler.method2, ...
            'dims', upscaler.dim2, 'nrelperm', upscaler.nrelperm, ...
            'pcow', false, 'verbose', false, 'cellinx', true);
        updata2  = upscaler2.upscale();
        data.krO = updata2.krO;
        data.krW = updata2.krW;
        
        if uspcaler.saveStep1
            % Save data from the first upscaling step, if requested
            data.krStep1 = updata1;
            data.krStep1.method    = upscaler.method1;
            data.krStep1.dim       = upscaler.dim1;
            data.krStep1.partition = p1;
        end
        
        if upscaler.verbose
            t = toc(t);
            fprintf('  Rel.perm #2:  %6.3f sec.\n', t);
            t = toc(t_relperm);
            fprintf('  Rel.perm tot: %6.3f sec.\n', t);
        end
        
        
        %------------------------------------------------------------------
        % Capillary pressure upscaling
        %------------------------------------------------------------------
        if upscaler.pcow
            t = tic;
            data = upPcOW(block, data, 'npointsmax', upscaler.npcow);
            if upscaler.verbose
                t = toc(t);
                fprintf('  Cap.pres:     %6.3f sec.\n', t);
            end
        end
        
    end
    
end

end


