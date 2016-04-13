function f = assignSURFADS(f, surfads, reg)
   surfads = extendTab(surfads);

   regmap = @(c, varargin) ...
      getRegMap(c, reg.SATNUM, reg.SATINX, varargin{:});

    f.surfads = @(c, varargin) interpReg(surfads, c, regmap(c, varargin{:}));
end
