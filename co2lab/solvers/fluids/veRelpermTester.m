function [s, pc, kr, SH, krH, s_max, fval] = veRelpermTester(hs, p, fluid, H, varargin)
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
   opt = struct('samples', 100, 'hs_max', [], 'kscale', [], 'const_res', false);
   opt = merge_options(opt, varargin{:});

   if(isfield(fluid, 'is_kscaled'))
      if(fluid.is_kscaled)
         invPc3D = @(p) fluid.invPc3D(p, opt.kscale);
      else
         invPc3D = @(p) fluid.invPc3D(p);
      end
   else
      invPc3D = @(p) fluid.invPc3D(p);
   end

   h = linspace(0, hs, opt.samples)';
   dh = h(2) - h(1);
   drho = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)) * norm(gravity);

   h_org = h(end) - h;
   s_h = invPc3D(h .* drho);
   sg_h = 1 - s_h;
   kr_h = fluid.kr3D(sg_h) .* (h_org<H);
   pc = hs .* drho;
   SH = (sum(sg_h .* (h_org<H)) - sg_h(end) .* (h_org(end)<H) / 2) .* dh;
   krH = (sum(kr_h) - kr_h(end) / 2) .* dh;
   if(~isempty(opt.hs_max))
      hs_max    = opt.hs_max;
      h_max     = linspace(0, hs_max, opt.samples)';
      h_max_org = h_max(end) - h_max;
      dh_max    = h_max(2) - h_max(1);
      s_hmax    = invPc3D(h_max .* drho);
      sg_hmax   = (1 - s_hmax) .* (h_max_org<H);
      SH_max    = (sum(sg_hmax) - sg_hmax(end) / 2) .* dh_max;
      s_max     = SH_max ./ H;

      if(~opt.const_res)
         SH = ((SH_max - SH) / (1 - fluid.res_water)) * fluid.res_gas + SH;
      else
         SH = SH + (min(opt.hs_max, H) - min(hs, H)) * fluid.res_gas;
      end

      s = SH ./ H;
   else
      s = SH ./ H;
   end
   kr = krH ./ H;
    fval = struct('h'         , h                 , ...
                  'h_org'     , h_org             , ...
                  's_h'       , s_h               , ...
                  'sg_h'      , sg_h .* (h_org<H) , ...
                  'kr_h'      , kr_h              , ...
                  's_hmax'    , s_hmax            , ...
                  'h_max'     , h_max             , ...
                  'h_max_org' , h_max_org         , ...
                  'sg_hmax'   , sg_hmax);
end