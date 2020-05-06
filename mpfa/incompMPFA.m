function state = incompMPFA(state, G, T, fluid, varargin)
    
    opt = struct('useTensorAssembly', false,...
                 'mpfastruct', []); 
    [opt, extra] = merge_options(opt, varargin{:});
    
    if ~opt.useTensorAssembly
        % use legacy implementation
        state = incompMPFAlegacy(state, G, T, fluid, varargin{:})
    else
        % use tensor assembly based implementation
        assert(~isempty(opt.mpfastruct), 'option mpfastruct is required for tensor assembly');
        mpfastruct = opt.mpfastruct;
        opt = struct('bc', [], ...
                     'W', []);
        [opt, extra] = merge_options(opt, extra{:});
        if ~isempty(opt.bc)
            bc = opt.bc;
            state = incompMPFAbc(G, mpfastruct, bc, extra{:});
        elseif ~isempty(opt.W)
            W = opt.W;
            state = incompMPFATensorAssembly(G, mpfastruct, W, extra{:});        
        else
            error('No stimulation input, either bc or W, given');
        end
    end
    


end

