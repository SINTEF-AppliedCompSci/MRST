function fsg = free_sg(sg, sGmax, rw, rn)
% Determine the mobile part of present saturation, for a hysteretic model
% with (potentially) a capillary fringe.  Homogeneous rock properties are 
% assumed. 
% 
% SYNOPSIS:
%   function fsg = free_sg(sg, sGmax, opt)
%
% DESCRIPTION:
% This function determine the amount of present saturation that is 
% effectively mobile.  It works for sharp-interface models, as well as 
% capillary fringe models that uses endpoint scaling.  'fsg' represent
% the CO2 that will be able to evacuate, e.g. excluding anything already
% residually trapped, or that will be left as residuall yrapped.
%
% As an example, for sharp-interface models (no capillary fringe), 
% the free saturation can be thought of the part of the CO2 saturation 
% that is within the mobile plume (as opposed to the residual saturation
% in the inebibided zone), and that is not part of what will be left behind 
% as residual trapping after imbibition.
%    
% For the sharp-interface model, the formula can be thought of as
% representing the simple transformation: 
%       s * H = h * (1 - sr(2)) + (h_max - h) * sr(1) 
%   s_max * H = h_max * (1 - sr(2))
% 
% For the capillary fringe model with endpoint scaling, the formula can be
% derived from the expression:
% s_free = s - (C/(1-C)) (s_max - s),  where C = rn / (1 - rw)
%
% Refer to formula 17 of the paper "Fully-implicit simulation of 
% vertical-equilibrium models with hysteresis and capillary fringe" 
% (Nilsen et al, 2016, DOI 10.1007/s10596-016-9547-y) for a derivation of the
% latter formula.
%   
% PARAMETERS:
%   sg    - present saturation
%   sGmax - historically maximum saturation
%   rw    - residual water saturation ('wetting saturation')
%   rn    - residual gas saturation   ('nonwetting saturation')
% RETURNS:
%   fsg - the free part of the present saturation (fsg <= sg)
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    %% Computing free part of current saturation
    % NB: in the h-formulation, this would equal (1-sw) h / H. 

    fsg = ((1 - rw) * sg - (sGmax * rn)) ./ (1 - rw - rn);

    %% Ensuring exact values at boundaries
    ineb = (sg >= sGmax);  % cells with no residual trapping (drainage zone)
    fsg = ifcond(sg, fsg, ineb);
    fsg = max(fsg, 0.0);
end 
