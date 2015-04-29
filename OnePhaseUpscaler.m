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
       
    function [data, report] = upscaleBlock(upscaler, block)
        
        data.dims = upscaler.dims;
        
        % Absolute permeability
        [data.perm, report] = upAbsPerm(block, 'dims', upscaler.dims, ...
            'dp', upscaler.dp);
        
        if upscaler.verbose
            fprintf('  Abs.perm:     % 2.3fs\n', report.time);
        end
        
        % Porosity
        data.poro = sum(block.pv) / sum(block.G.cells.volumes);
        
    end
    
end

end

