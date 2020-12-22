function fluid = include_BO_form(fluid, shortname, ref_val)

   % Add reference value
   fluid.(['rho', shortname, 'S']) = ref_val;

   rhofun = fluid.(['rho', shortname]);

   if nargin(rhofun) > 1
      % density is considered function of P and T
      bfun = @(p, t, varargin) rhofun(p, t) ./ ref_val;
   else
      % density is considered a function of P only
      bfun = @(p, varargin) rhofun(p) ./ ref_val;
   end

   % Add fluid formation factor
   fluid.(['b', shortname]) = bfun;

end