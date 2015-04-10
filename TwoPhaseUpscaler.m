classdef TwoPhaseUpscaler < OnePhaseUpscaler
%Two phase upscaling

properties
    method
    nrelperm
    pcow
    npcow
end

methods
    
    function upscaler = TwoPhaseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid);
        
        upscaler.method   = [];
        upscaler.nrelperm = 20;
        upscaler.pcow     = true; % Upscale pcow or not
        upscaler.npcow    = 100;
        upscaler = merge_options(upscaler, varargin{:});
    end
       
    function data = upscaleBlock(upscaler, block)
        
        % Perform one phase upscaling first
        data = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        
        % Capillary pressure upscaling
        if upscaler.pcow
            t = tic;
            data = upPcOW(block, data, 'npointsmax', upscaler.npcow);
            if upscaler.verbose
                t = toc(t);
                fprintf('  Cap.pres:     %6.3f sec.\n', t);
            end
        end
        
        % Relative permeability upscaling
        t = tic;
        data = upRelPerm(block, data, upscaler.method, ...
            'nvalues', upscaler.nrelperm, 'dims', upscaler.dims, ...
            'dp', upscaler.dp);
        if upscaler.verbose
            t = toc(t);
            fprintf('  Rel.perm:     %6.3f sec.\n', t);
        end
        
    end
    
end

end


