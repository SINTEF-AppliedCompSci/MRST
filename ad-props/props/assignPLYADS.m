function f = assignPLYADS(f, plyads, reg)
   plyads = extendTab(plyads);

   regmap = @(c, varargin) ...
      getRegMap(c, reg.SATNUM, reg.SATINX, varargin{:});

    f.ads = @(c, varargin) interpReg(plyads, c, regmap(c, varargin{:}));
end
