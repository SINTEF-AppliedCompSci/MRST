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
        switch upscaler.OnePhaseMethod
            case 'pressure'
                [data.perm, report] = upAbsPerm(block, ...
                    'dims', upscaler.dims, 'dp', upscaler.dp);
            case {'arithmetic', 'harmonic', 'geometric', ...
                    'harmonic-arithmetic'}
                [data.perm, report] = upAbsPermAvg(block, ...
                    'dims', upscaler.dims, ...
                    'method', upscaler.OnePhaseMethod);
            otherwise
                error('Unknown one-phase upscaling method.');
        end
        
        if upscaler.verbose
            fprintf('  Abs.perm:     % 2.3fs\n', report.time);
        end
        
        % Porosity
        data.poro = sum(block.pv) / sum(block.G.cells.volumes);
        
    end
    
end

end

