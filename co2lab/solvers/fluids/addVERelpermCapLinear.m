function fluid = addVERelpermCapLinear(fluid, cap_scale, varargin)
% VE relperm with linear capillary pressure

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
   opt = struct('res_water', 0, 'res_gas', 0, 'beta', 1, 'H', [], 'kr_pressure', false); 
   opt = merge_options(opt, varargin{:}); 
   
   % should also include endpoint scaling
   
   end_scale = (1 - opt.res_water).^opt.beta; 
   if(~opt.kr_pressure)
      fake_pressure = 200 * barsa; % @@ Make the value of this constant optional?
      fluid.krG = @(sg, p, varargin) end_scale .* krG(sg, fake_pressure, fluid, ...
                                                      cap_scale, opt, varargin{:}); 
      fluid.krW = @(sw, p, varargin) krW(sw, opt, varargin{:}); 
      fluid.pcWG = @(sg, p, varargin) pcWG(sg, p, fluid, cap_scale, opt, varargin{:}); 
      % fluid.S3D = @(SVE, samples, H) S3D(SVE, fake_pressure, samples, H, fluid, ...
      %                                    cap_scale, opt); 
      
   else
      fluid.krG = @(sg, p, varargin) end_scale .* krG(sg, p, fluid, cap_scale, ...
                                                      opt, varargin{:}); 
      fluid.krW = @(sw, p, varargin) krW(sw, opt, varargin{:}); 
      fluid.pcWG = @(sg, p, varargin) pcWG(sg, p, fluid, cap_scale, opt, varargin{:}); 
   end
   fluid.res_water = opt.res_water; 
   fluid.res_gas = opt.res_gas; 
   fluid.invPc3D = @(p) invPc(p, cap_scale, opt); 
   fluid.kr3D = @(s) end_scale .* (s ./ (1 - opt.res_water)).^opt.beta; 
end

function s = invPc(p, cap_scale, opt)
   s = (p / cap_scale) .* (1 - opt.res_water); 
   s(s > (1 - opt.res_water)) = 1 - opt.res_water; 
   s(s < 0) = 0; 
   s = 1 - s; 
end

function kr = krG(sg, p, fluid, cap_scale, opt, varargin)
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   if(~isempty(loc_opt.sGmax))
      sg_free = free_sg(sg, loc_opt.sGmax, opt); 
      h = invS(sg_free, p, fluid, cap_scale, opt); 
      kr = S_beta(h, p, opt.beta, fluid, cap_scale, opt); 
   else
      h = invS(sg, p, fluid, cap_scale, opt); 
      kr = S_beta(h, p, opt.beta, fluid, cap_scale, opt); 
   end
end

function kr = krW(sw, opt, varargin)
% beta = 1; 
% fo now linear
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   if(~isempty(loc_opt.sGmax))
      sg = 1 - sw; 
      ineb = (sg)>loc_opt.sGmax; 
      sg_res = (loc_opt.sGmax - sg); 
      % assert(all(sg_res >= 0)); 
      sw_free = 1 - loc_opt.sGmax; 
      kr = sw_free + (1 - opt.res_gas) * sg_res; 
      % this to avoid errors in ADI derivative
      kr(~ineb) = sw(~ineb); 
      kr(kr<0) = 0.0 * kr(kr<0); 
      assert(all(kr >= 0)); 
   else
      kr = sw;
    end
end

function pc = pcWG(sg, p, fluid, cap_scale, opt, varargin)

   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 

   if(~isempty(loc_opt.sGmax))
      sg_free = free_sg(sg, loc_opt.sGmax, opt); 
      h = invS(sg_free, p, fluid, cap_scale, opt); 
      assert(all(sg_free >= 0))
      pc = norm(gravity) * (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)) .* h; 
   else
      h = invS(sg, p, fluid, cap_scale, opt); 
      pc = norm(gravity) * (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)) .*h;
   end
end

function S_out = S_beta(h, p, beta, fluid, cap_scale, opt)
% integral of power of Saturation
% @@ There seems to be differences between this formula and the one expressed
% in section 3.2 of paper 3.  Check!
   C = cap_scale; 
   drho = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)); 
   cap_norm = (C ./ (drho * norm(gravity))) ./ opt.H; 
   % V_cap = cap_norm.^(beta + 1) / (beta + 1); 
   % V_cap = cap_norm / (beta + 1); 
   h_norm = h ./ opt.H; 
   S_out = (h_norm - cap_norm) + cap_norm ./ (beta + 1); 
   ind = h_norm<cap_norm; 
   S_out(ind) = (h_norm(ind) ./ cap_norm(ind)).^(beta + 1) / (beta + 1) .* cap_norm(ind); 
   ind = h_norm>1; 
   if(any(ind))
      S_out(ind) = S_out(ind) - ((h_norm(ind) - 1) ./ cap_norm(ind)).^(beta +1)/(beta+1) .* cap_norm(ind);
    end
end

function h = invS(S, p, fluid, cap_scale, opt)
   
   S = S / (1 - opt.res_water); 
   C = cap_scale; 
   drho = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)); 
   cap_norm = (C ./ (drho * norm(gravity))) ./ opt.H; 
   vol_cap = (1 / 2) * cap_norm; 
   assert(all(cap_norm<1)); 
   vol_capb = vol_cap + 1 - cap_norm; 
   h = (2 * (S) ./ cap_norm).^(1 / 2) .* cap_norm; 
   
   % regularize square root function near zero
   h_tol = 1e-5; 
   S_tol = 0.5 * (h_tol ./ value(cap_norm)).^2 .* value(cap_norm); 
   h_tmp = h_tol * S ./ S_tol; 

   ind = h<h_tol; 
   h(ind) = h_tmp(ind); 
   
   ind = (S>vol_cap) & (S<vol_capb); 
   if(any(ind))
      h(ind) = (S(ind) - vol_cap(ind)) + cap_norm(ind); 
   end

   ind = (S >= vol_capb); 
   if(any(ind))
      a = -1 ./ (2 * cap_norm(ind)); b = 1; c = -S(ind) + (vol_cap(ind) - cap_norm(ind) + 1); 
      res = (b).^2 - 4 * a .* c; 
      res(res<0) = 0; 
      h(ind) = ( - b + res.^(1 / 2)) ./ (2 * a) + 1; 

   end
   h = h .* opt.H; 

 end
 
 
% function [s, z_out] = S3D(SVE,p, samples, H, fluid, cap_scale, opt)
%     % utility function for plotting saturation frofile and definition of
%     % the saturation profile the other functions are obtained from    
%     C=cap_scale;  
%     drho  =  (fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p));
%     cap_norm=(C./(drho*norm(gravity)))./H;
%     %vol_cap=cap_norm.^2/2;
%     vol_cap=1/2*cap_norm;
%     vol_capb=vol_cap+1-cap_norm;
%     h=invS(S,p,fluid,opt);
%     z=linspace(0,h,samples);
%     if(SVE<vol_cap)
%        s=z*(h./(opt.H*cap_norm));
%     else (SVE<vol_capb)
%        s=1; 
%        s(z<cap_norm*opt.H)=z./(cap_norm*opt.H); 
%     end
%     z_out=h-z;
%     s=s(end:-1:1);
% end
 