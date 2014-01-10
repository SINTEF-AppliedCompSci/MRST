function f = assignROCKTAB(f, rocktab, reg)
   if ~iscell(rocktab), rocktab = { rocktab }; end

   Tpv = cellfun(@(x) x(:, [1, 2]), rocktab, 'UniformOutput', false);
   Ttr = cellfun(@(x) x(:, [1, 3]), rocktab, 'UniformOutput', false);

   regmap = @(p, varargin) ...
      getRegMap(p, reg.ROCKNUM, reg.ROCKINX, varargin{:});

   f.pvMultR   = @(p, varargin) interpReg(Tpv, p, regmap(p, varargin{:}));
   f.tranMultR = @(p, varargin) interpReg(Ttr, p, regmap(p, varargin{:}));
end
