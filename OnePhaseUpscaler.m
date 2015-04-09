classdef OnePhaseUpscaler < Upscaler
%One phase upscaling

properties
    dp
    dims
end

methods
    
    function upscaler = OnePhaseUpscaler(G, rock, varargin)
        upscaler = upscaler@Upscaler(G, rock);
        
        upscaler.dp   = 1*barsa;
        upscaler.dims = 1:3;
        upscaler = merge_options(upscaler, varargin{:});
    end
       
    function data = upscaleBlock(upscaler, block)
        data.K = upAbsPerm(block, 'dims', upscaler.dims, ...
            'dp', upscaler.dp);
    end
    
end

end

