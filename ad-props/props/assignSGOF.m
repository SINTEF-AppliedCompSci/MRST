function f = assignSGOF(f, sgof, reg)
   cfun  = @(f) cellfun(f, sgof, 'UniformOutput', false);

   % Compute tables (static data)
   Tkrg  = extendTab( cfun(@(x) x(:, [1, 2])) );
   Tkro  = extendTab( cfun(@(x) x(:, [1, 3])) );
   Tpcog = extendTab( cfun(@(x) x(:, [1, 4])) );

   % Region mapping
   regmap = @(sw, varargin) ...
      getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});

   if ~ isfield(f, 'sWcon'),
      sgas = @(so) 1 - so;
   else
      sgas = @(so) 1 - so - f.sWcon;
   end

   % Region interpolator
   ireg = @(T, sg, varargin) interpReg(T, sg, regmap(sg, varargin{:}));

   f.krG  = @(sg, varargin) ireg(Tkrg , sg      , varargin{:});
   f.krOG = @(so, varargin) ireg(Tkro , sgas(so), varargin{:});
   f.pcOG = @(sg, varargin) ireg(Tpcog, sg      , varargin{:});
end
