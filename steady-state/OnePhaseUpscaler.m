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
