classdef TwoPhaseTwoStepUpscaler < OnePhaseUpscaler
    %Two phase upscaling

properties
    RelpermMethod1
    RelpermMethod2
    RelpermAbsMethod1
    RelpermAbsMethod2
    dim1
    dim2
    nrelperm
    pcow
    npcow
    gravity
    
    saveStep1
    twostepz
end

methods
    
    function upscaler = TwoPhaseTwoStepUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid);
        
        upscaler.RelpermMethod1   = [];
        upscaler.RelpermMethod2   = [];
        upscaler.RelpermAbsMethod1 = 'pressure';
        upscaler.RelpermAbsMethod2 = 'pressure';
        upscaler.dim1      = [];
        upscaler.dim2      = [];
        upscaler.nrelperm  = 20;
        upscaler.pcow      = true; % Upscale pcow or not
        upscaler.npcow     = 100;
        upscaler.gravity   = false;
        
        upscaler.saveStep1 = false; % save data from first upscaling
        upscaler.twostepz  = false;
        
        upscaler = merge_options(upscaler, varargin{:});
        
        % The absolute permeability will be evaluated in the second
        % dimension chosen.
        % upscaler.dims = upscaler.dim2;
        
    end
       
    function [data, report] = upscaleBlock(upscaler, block)
        
        % Because the slizes will contain only a single cell in a
        % direction, periodic boundaries cannot be applied without using
        % some tricks
        assert(~block.periodic, ['The two-stage upscaling does not '...
            'support periodic boundary']);
        
        wantReport = nargout > 1;
        startTime  = tic;
        
        % Separate special upscaling in the z-direction is implemented.
        upscalez = upscaler.twostepz;
        
        assert(~isempty(upscaler.RelpermMethod1) && ...
            ~isempty(upscaler.RelpermMethod2), 'Methods must be given.');
        assert(~isempty(upscaler.dim1) && ...
            ~isempty(upscaler.dim2), 'Dimensions must be given.');
        
        % Perform one phase upscaling first
        [data, rep] = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        if wantReport
            report.onephase = rep;
        end
        
        %------------------------------------------------------------------
        % Relative permeability upscaling, step 1
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
            'partition', p1, 'OnePhaseMethod', upscaler.OnePhaseMethod, ...
            'RelpermMethod', upscaler.RelpermMethod1, ...
            'RelpermAbsMethod', upscaler.RelpermAbsMethod1, ...
            'dims', dims1, 'nrelperm', upscaler.nrelperm, ...
            'pcow', true, 'verbose', false, 'deck', block.deck);
        if wantReport
            [updata1, ~, rep] = upscaler1.upscale();
            report.relperm.method = 'twostep';
            report.relperm.step1  = rep;
        else
            updata1 = upscaler1.upscale();
        end
        relperm1Time = toc(t);
        
        if upscaler.verbose
            fprintf('  Rel.perm #1:  %6.3fs\n', relperm1Time);
        end
        
        
        %------------------------------------------------------------------
        % Relative permeability upscaling, step 2
        % (Upscale in direction 2, using method 2)
        %------------------------------------------------------------------
        
        % Direction 2 is typically x or y
        
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
            'partition', p2, 'OnePhaseMethod', upscaler.OnePhaseMethod, ...
            'RelpermMethod', upscaler.RelpermMethod2, ...
            'RelpermAbsMethod', upscaler.RelpermAbsMethod2, ...
            'dims', upscaler.dim2, 'nrelperm', upscaler.nrelperm, ...
            'pcow', false, 'verbose', false);
        if upscalez
            upscaler2.savesat = true; % save sat. distributions
            [updata2, ~, rep] = upscaler2.upscale();
            satdist = rep.blocks{1}.relperm.satdist;
            rep.blocks{1}.relperm = rmfield(rep.blocks{1}.relperm, ...
                'satdist');
        elseif wantReport
            [updata2, ~, rep] = upscaler2.upscale();
        else
            updata2 = upscaler2.upscale();
        end
        if wantReport
            report.relperm.step2 = rep;
        end
        krO = updata2.krO;
        krW = updata2.krW;
        
        % Upscale relperm in z-direction in a special manner if requested.
        if upscalez
            sW   = krO{1}(:,1);
            nsat = numel(sW);
            krWz = nan(nsat, 1);
            krOz = nan(nsat, 1);
            
            % interp function
            int  = @(T,x) interp1(T(:,1), T(:,2), x, 'linear', 'extrap');
            
            % Loop over saturation distributions found in the first step of
            % the upscaling. That is, distributions on the intermediate
            % coarse grid.
            for is=1:numel(satdist)
                dist = satdist{is};
                
                % Find the relperm value in each cell of the inter. grid
                krWi = nan(nsat, 1);
                krOi = nan(nsat, 1);
                for ic=1:numel(dist) % loop over cells
                    krWi(ic) = int(updata1(ic).krW{2}, dist(ic));
                    krOi(ic) = int(updata1(ic).krO{2}, dist(ic));
                end
                
                % The upscaled relperm value is taken as the arithmetic
                % mean of the block relperms
                krWz(is) = mean(krWi);
                krOz(is) = mean(krOi);
            end
            
            % Ensure endpoints are at zero
%             krWz(1)   = 0;
%             krOz(end) = 0;
            
            % Using x-direction also in y-direcion (or vica versa)
            krW = {krW{1}, krW{1}, [sW krWz]};
            krO = {krO{1}, krO{1}, [sW krOz]};
        end
        
        % Save to structure
        data.krO = krO;
        data.krW = krW;
        relperm2Time = toc(t);
        
        if upscaler.saveStep1
            % Save data from the first upscaling step, if requested
            data.krStep1 = updata1;
            data.krStep1.method    = upscaler.RelpermMethod1;
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
    
    
    function block = createBlock(upscaler, cells)
    % Create grid, rock and fluid for the sub block represented by the
    % given cells.
        
        block = createBlock@OnePhaseUpscaler(upscaler, cells);
        
        if ~isempty(block.deck)
            
            % To use capillary upscaling, we need a function pcOWInv
            block.fluid = addPcOWInvADIFluid(block.fluid, block.deck);
            
            % To use viscous upscaling, we need a function fracFlowInv
            block.fluid = addFracFlowInvADIFluid(block.fluid, block.deck);
            
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
