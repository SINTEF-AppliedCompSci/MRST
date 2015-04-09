classdef OnePhaseUpscaler < Upscaler
%One phase upscaling

properties
    dp
    boundaryFaces
end

methods
    
    function upscaler = OnePhaseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@Upscaler(G, rock, fluid);
    end
       
    function data = upscaleBlock(upscaler, G, rock, fluid)
        %error('Needs to be overridden');
        fprintf('ONE PHASE Single block upscaling of %d cells.\n', numel(G.cells.num));
        data.val = 1;
    end
    
    function [V, states] = findFluxes(G, rock)

        dims  = opt.dims;
        ndims = length(dims);

        % Initial state
        state0 = initResSol(G, 100*barsa, 1);

        % Apply pressure drop in each dimension
        for i = 1:ndims
            
            % Set boundary conditions
            bc = addBC([], upscaler.boundaryFaces{dims(i),1}, 'pressure', upscaler.dp);
            bc = addBC(bc, upscaler.boundaryFaces{dims(i),2}, 'pressure', 0);

            % Solve
            warning('off','mrst:periodic_bc');
            state1 = opt.psolver(state0, rock, bc);
            warning('on','mrst:periodic_bc');

            faces = opt.bcFaces{dims(i), 2};
            sign  = ones(numel(faces), 1);
            sign(G.faces.neighbors(faces,1)==0) = -1;
            V(i) = sum(state1.flux(faces, 1).*sign) / opt.areas(i);
        end
    
    end
    
end

end

