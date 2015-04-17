function f = assignPVDO(f, pvdo, reg)
   cfun    = @(f) cellfun(f, pvdo, 'UniformOutput', false);

   % Compute tables (static data)
   TBO     = cfun(@(x) x(:, [1, 2]));
   TbO     = cfun(@(x) [x(:,1), 1 ./ x(:,2)]);
   TmuO    = cfun(@(x) x(:, [1, 3]));
   TBOxmuO = cfun(@(x) [x(:,1), prod(x(:, [2, 3]), 2)]);

   % Region mapping
   regmap = @(po, varargin) ...
      getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});

   % Region interpolator
   ireg = @(T, po, varargin) interpReg(T, po, regmap(po, varargin{:}));

   f.BO     = @(po, varargin) ireg(TBO,     po, varargin{:});
   f.bO     = @(po, varargin) ireg(TbO,     po, varargin{:});
   f.muO    = @(po, varargin) ireg(TmuO,    po, varargin{:});
   f.BoxmuO = @(po, varargin) ireg(TBOxmuO, po, varargin{:});
end
