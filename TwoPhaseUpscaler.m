classdef TwoPhaseUpscaler < OnePhaseUpscaler
%Two phase upscaling

properties
    method
    nrelperm
    pcow
    npcow
    
    relpermdims  % Dimensions to upscale relperm in
    savesat      % Save saturation distributions
    
end

methods
    
    function upscaler = TwoPhaseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid, ...
            varargin{:});
        
        upscaler.method   = [];
        upscaler.nrelperm = 20;
        upscaler.pcow     = true; % Upscale pcow or not
        upscaler.npcow    = 100;
        upscaler.relpermdims = upscaler.dims;
        upscaler.savesat  = false; % Save saturation distributions
        
        upscaler = merge_options(upscaler, varargin{:});
    end
       
    function [data, report] = upscaleBlock(upscaler, block)
        
        wantReport  = nargout > 1;
        startTime   = tic;
        
        % Perform one phase upscaling first
        [data, rep] = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        if wantReport
            report.onephase = rep;
        end
        
        data.relpermdims = upscaler.relpermdims;
        
        % Capillary pressure upscaling
        % If using pore volume relperm upscaling, we also upscale the pcOW
        % in the relperm function, so don't do it here.
        if upscaler.pcow && ~strcmpi(upscaler.method, 'porevolume')
            
            if ~isempty(regexpi(upscaler.method, '_grav$'));
                % We upscale with gravity
                [data, rep] = upPcOW(block, data, ...
                    'npointsmax', upscaler.npcow, 'gravity', 'centroid');
                % We also need the bottom pcow curve for relperm upscaling
                [d, r] = upPcOW(block, [], ...
                    'npointsmax', upscaler.npcow, 'gravity', 'bottom');
                data.pcOW_bot = d.pcOW;
                if wantReport
                    report.pcow_bot = r;
                end
            else
                % Default: no gravity
                [data, rep] = upPcOW(block, data, ...
                    'npointsmax', upscaler.npcow);
            end
            if wantReport
                report.pcow = rep;
            end
            if upscaler.verbose
                fprintf('  Cap.pres:     %6.3fs\n', rep.time);
            end
        end
        
        % Relative permeability upscaling
        if strcmpi(upscaler.method, 'porevolume')
            up = @() upRelPermPV(block, data, upscaler.method, ...
                'nvalues', upscaler.nrelperm);
        else
            up = @() upRelPerm(block, data, upscaler.method, ...
                'nvalues', upscaler.nrelperm, 'dp', upscaler.dp, ...
                'dims', upscaler.relpermdims, 'savesat', upscaler.savesat);
        end
        if wantReport
            [data, rep] = up();
            report.relperm = rep;
        else
            data = up();
        end
        if upscaler.verbose
            fprintf('  Rel.perm:     %6.3fs\n', rep.time);
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


