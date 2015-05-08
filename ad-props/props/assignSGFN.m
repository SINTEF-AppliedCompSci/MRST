function f = assignSGFN(f, sgfn, reg)
   cfun = @(f) cellfun(f, sgfn, 'UniformOutput', false);

   % Compute tables (static data)
   Tkrg  = extendTab( cfun(@(x) x(:, [1, 2])) );
   Tpcog = extendTab( cfun(@(x) x(:, [1, 3])) );

   % Region mapping
   regmap = @(sg, varargin) ...
      getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});

   % Region interpolator
   ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));

   f.krG  = @(sg, varargin) ireg(Tkrg , sg, varargin{:});
   f.pcOG = @(sg, varargin) ireg(Tpcog, sg, varargin{:});
end
