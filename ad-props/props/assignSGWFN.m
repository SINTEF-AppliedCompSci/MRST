function f = assignSGWFN(f, sgfn, reg)
   cfun = @(f) cellfun(f, sgfn, 'UniformOutput', false);

   % Compute tables (static data)
   Tkrg  = extendTab( cfun(@(x) x(:, [1, 2])) );
   Tkrw  = extendTab( cfun(@(x) flipud([1 - x(:, 1), x(:, 3)])));
   Tpcwg = extendTab( cfun(@(x) x(:, [1, 4])) );

   % Region mapping
   regmap = @(sg, varargin) ...
      getRegMap(sg, reg.SATNUM, reg.SATINX, varargin{:});

   % Region interpolator
   ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));

   f.krG  = @(sg, varargin) ireg(Tkrg , sg, varargin{:});
   f.krW  = @(sw, varargin) ireg(Tkrw, sw, varargin{:});
   f.pcWG = @(sg, varargin) ireg(Tpcwg, sg, varargin{:});
end