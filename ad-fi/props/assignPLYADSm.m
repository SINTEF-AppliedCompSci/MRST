function f = assignPLYADSm(f, plyvisc, reg)
   regmap = @(c, varargin) ...
      getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});

   f.muWMult = @(c, varargin) ...
      interpReg(plyvisc, c, regmap(c, varargin{:}));
end
