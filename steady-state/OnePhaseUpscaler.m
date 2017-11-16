classdef OnePhaseUpscaler < Upscaler
    %One phase upscaling

properties
    dp
    dims
    OnePhaseMethod
end

methods
    
    function upscaler = OnePhaseUpscaler(G, rock, varargin)
        upscaler = upscaler@Upscaler(G, rock, varargin{:});
        
        upscaler.dp   = 1*barsa;
        upscaler.dims = 1:3;
        upscaler.OnePhaseMethod = 'pressure';
        upscaler = merge_options(upscaler, varargin{:});
    end
       
    function [data, report] = upscaleBlock(upscaler, block)
        
        data.dims = upscaler.dims;

        % Absolute permeability
        [data, report] = upAbsPerm(block, data, ...
            'method', upscaler.OnePhaseMethod, ...
            'dims', upscaler.dims, ...
            'dp', upscaler.dp);
        
        if upscaler.verbose
            fprintf('  Abs.perm:     % 2.3fs\n', report.time);
        end
        
        % Porosity
        data = upPoro(block, data);
        
    end
    
end

end

