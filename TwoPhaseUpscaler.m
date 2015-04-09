classdef TwoPhaseUpscaler < OnePhaseUpscaler
%Two phase upscaling

properties
    method
    nvalues
end

methods
    
    function upscaler = TwoPhaseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid);
        
        upscaler.method  = [];
        upscaler.nvalues = 20;
        upscaler = merge_options(upscaler, varargin{:});
    end
       
    function data = upscaleBlock(upscaler, block)
        
        % Perform one phase upscaling first
        data = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        
        % Two phsase upscaling
        data = upRelPerm(block, data, upscaler.method, ...
            'nvalues', upscaler.nvalues, 'dims', upscaler.dims, ...
            'dp', upscaler.dp);
        
    end
    
end

end


