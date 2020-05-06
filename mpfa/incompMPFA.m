function state = incompMPFA(state, G, mpfaT, fluid, varargin)
% Solve incompressible problem with mpfa transmissibilities.
%
% Two versions available : 'legacy' (default) and 'tensor Assembly'.
%
% This legacy version is faster. It is limited to a mesh with grid cells where
% corners have the same number of faces as the spatial dimension (this is always
% the case in 2D but not in 3D). The tensor assembly version
% (computeMultiPointTransTensorAssembly) can handle the other cases but is
% slower (the implementation will be optimized in the future to run faster).
    opt = struct('useTensorAssembly', false,...
                 'mpfastruct', []); 
    [opt, extra] = merge_options(opt, varargin{:});
    
    if ~opt.useTensorAssembly
        
        % use legacy implementation
        state = incompMPFAlegacy(state, G, mpfaT, fluid, varargin{:})
    
    else
        
        % use tensor assembly based implementation
        state = incompMPFATensorAssembly(G, mpfaT, extra{:});        

    end
    


end

