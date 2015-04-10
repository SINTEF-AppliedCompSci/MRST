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
        
        % Absolute permeability
        t = tic;        
        data.perm = upAbsPerm(block, 'dims', upscaler.dims, ...
            'dp', upscaler.dp);
        if upscaler.verbose
            t = toc(t);
            fprintf('  Abs.perm:     % 2.3fs\n', t);
        end
        
        % Porosity
        data.poro = sum(block.pv) / sum(block.G.cells.volumes);
        
    end
    
end

end

