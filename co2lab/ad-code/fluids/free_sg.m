function fsg = free_sg(sg, sGmax, opt)
% Determine the mobile part of present saturation.
% 
% SYNOPSIS:
%   function fsg = free_sg(sg, sGmax, opt)
%
% DESCRIPTION:
% Determine the mobile part of present saturation, given the present saturation
% and the historically maximum saturation, assuming a sharp interface.  (If the
% present saturation is lower than the historical maximum, it suggest that some
% of it is residually trapped, leading to a lower mobile saturation).  The
% formula is based on the simple transformation: 
%     s * H = h * (1 - sr(2)) + (h_max - h) * sr(1)
% s_max * H = h_max * (1 - sr(2))
% 
% PARAMETERS:
%   sg    - present saturation
%   sGmax - historically maximum saturation
%   opt   - structure expected to contain the following fields:
%           * opt.res_gas : residual gas saturation
%           * opt.res_water : residual oil saturation
% RETURNS:
%   fsg - the free part of the present saturation (fsg <= sg)

    %% aliases for better readability
    rw = opt.res_water; % 'wetting phase'
    rn = opt.res_gas; % 'non-wetting phase'
    
    %% Computing free part of current saturation
    % NB: in the h-formulation, this would equal (1-sw) h / H. 

    fsg = ((1 - rw) * sg - (sGmax * rn)) ./ (1 - rw - rn);

    %% Ensuring exact values at boundaries
    ineb = (sg >= sGmax);  % cells with no residual trapping (drainage zone)
    fsg = ifcond(sg, fsg, ineb);
    fsg = max(fsg, 0.0);
end 
