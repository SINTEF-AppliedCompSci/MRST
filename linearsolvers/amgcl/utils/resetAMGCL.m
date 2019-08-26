function resetAMGCL(varargin)
    amg_opt = getAMGCLMexStruct(varargin{:});
    amgcl_matlab(sparse([], [], []), [], amg_opt, nan, nan, 1000);
end