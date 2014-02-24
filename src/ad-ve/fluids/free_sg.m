function fsg = free_sg(sg, sGmax, opt)
% Determine the mobile part of present saturation, given the present
% saturation and the historically maximum saturation, assuming a sharp
% interface.   (If the present saturation is lower than the historical
% maximum, it suggest that some of it is residually trapped, leading to a
% lower mobile saturation).  The formula is based on the simple
% transformation:
% s*H=h*(1-sr(2))+(h_max -h)*sr(1)
% s_max*H = h_max*(1-sr(2))
% 
% SYNOPSIS:
%   function fsg = free_sg(sg, sGmax, opt)
%
% PARAMETERS:
%   sg    - present saturation
%   sGmax - historically maximum saturation
%   opt   - structure expected to contain the following fields:
%           * opt.res_gas : residual gas saturation
%           * opt.res_oil : residual oil saturation
% RETURNS:
%   fsg - the free part of the present saturation (fsg <= sg)
%
        ineb=sg>=sGmax;
        %sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;                
        fsg=((1-opt.res_oil)*sg-(sGmax*opt.res_gas))./(1-opt.res_gas-opt.res_oil);
        %fsg=(sg-(sGmax*opt.res_gas))./(1-opt.res_gas);
        %fsg(sg>=sGmax)=sg(sg>=sGmax);
        %assert(all(fsg>=-sqrt(eps)));
        fsg(ineb)=sg(ineb);
        fsg(fsg<0)=0.0*fsg(fsg<0);
end 
    
    % this transformation is based on the simple transormation
   % s*H=h*(1-sr(2))+(h_max -h)*sr(1)
   % s_max*H = h_max*(1-sr(2))
   