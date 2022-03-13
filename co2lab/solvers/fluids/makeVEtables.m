function tables = makeVEtables(varargin)
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

opt = struct('invPc3D'    , []    , ...
                'is_kscaled' , false , ...
                'kr3D'       , []    , ...
                'Gt'         , []    , ...
                'samples'    , 1000  , ...
                'drho'       , 400   , ...
                'kscale'     , []);
   opt = merge_options(opt, varargin{:});

   % make tables using h
   %H_max = 10 * max(opt.Gt.cells.H);
   H_max = 20 * max(opt.Gt.cells.H);
   
   h = linspace(0, H_max, opt.samples)';
   dh = h(2) - h(1);
   s_h = opt.invPc3D(h * opt.drho * norm(gravity));
   sg_h = 1 - s_h;
   kr_h = opt.kr3D(sg_h);
   
   %
   SH = (cumsum(sg_h) - sg_h / 2) * dh;
   SH_org = SH;
   %%
   ddh = h - mean(opt.Gt.cells.H);
   SH(ddh>0) = SH(ddh>0) - interpTable(h, SH_org, ddh(ddh>0));
   krH = (cumsum(kr_h) - kr_h / 2) * dh;
   krH_org = krH;
   krH(ddh>0) = krH(ddh>0) - interpTable(h, krH_org, ddh(ddh>0));

   if(opt.is_kscaled)
      p_max = max(H_max * opt.drho * norm(gravity) ./ opt.kscale);
   else
      p_max = H_max * opt.drho * norm(gravity);
   end
   
   p = linspace(0, p_max, opt.samples)';
   dp = p(2) - p(1);
   s_p = opt.invPc3D(p);
   sg_p = 1 - s_p;
   kr_p = opt.kr3D(sg_p);
   % simpson
   SP = (cumsum(sg_p) - sg_p / 2) * dp;
   krP = (cumsum(kr_p) - kr_p / 2) * dp;

   [p, SP, krP] = cleanTables(p, SP, krP);
   [h, SH, krH] = cleanTables(h, SH, krH);

   if(~opt.is_kscaled)
      tables = struct('p', p, 'SP', SP, 'krP', krP,...
                      'h', h, 'SH', SH, 'krH', krH,...
                      'is_kscaled', opt.is_kscaled,...
                      'H_max', H_max, 'drho', opt.drho);
   else
      tables = struct('p', p, 'SP', SP, 'krP', krP,...
                      'is_kscaled', opt.is_kscaled,...
                      'P_max', p_max, 'drho', opt.drho);
   end
   tables.invPc3D = @(p) opt.invPc3D(p);
   tables.kr3D = @(s)    opt.kr3D(s);
end

% ----------------------------------------------------------------------------


function [h, SH, krH] = cleanTables(h, SH, krH)
   [SH, ia] = unique(SH);
   h = h(ia);
   krH = krH(ia);
   assert(all(diff(SH)>0));
   ind = diff(SH)>1e-5;
   ind = [true; ind];
   h = h(ind);
   SH = SH(ind);
   krH = krH(ind);
end