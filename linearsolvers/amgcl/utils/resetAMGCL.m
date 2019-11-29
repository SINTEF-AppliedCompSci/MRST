function resetAMGCL(varargin)
    amg_opt = getAMGCLMexStruct(varargin{:});
    amgcl_matlab(sparse([], [], []), zeros(0, 1), amg_opt, nan, nan, 1000);
end