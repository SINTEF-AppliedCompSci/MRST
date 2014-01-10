function f = assignPVDG(f, pvdg, reg)
   cfun = @(f) cellfun(f, pvdg, 'UniformOutput', false);

   % Compute tables (static data)
   TBG  = cfun(@(x) x(:, [1, 2]));
   TbG  = cfun(@(x) [x(:,1), 1 ./ x(:,2)]);
   TmuG = cfun(@(x) x(:, [1, 3]));

   % Region mapping
   regmap = @(pg, varargin) ...
      getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});

   % Region interpolator
   ireg = @(T, pg, varargin) interpReg(T, pg, regmap(pg, varargin{:}));

   f.BG  = @(pg, varargin) ireg(TBG , pg, varargin{:});
   f.bG  = @(pg, varargin) ireg(TbG , pg, varargin{:});
   f.muG = @(pg, varargin) ireg(TmuG, pg, varargin{:});
end
