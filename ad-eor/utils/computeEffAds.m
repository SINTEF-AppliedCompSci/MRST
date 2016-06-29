function ads = computeEffAds(c, adsmax, fluid)
%
%
% SYNOPSIS:
%   function ads = computeEffAds(c, adsmax, fluid)
%
% DESCRIPTION: Compute the value of the adsorption, given the concentration
% and, if no desorption is allowed, a maximum adsorption value.
%
% PARAMETERS:
%   c      - Surfactant concentration
%   adsmax - Maximum adsorption value
%   fluid  - Fluid structure
%
% RETURNS:
%   ads - Adsorption
%
% EXAMPLE:
%
% SEE ALSO:
%
    ads = fluid.surfads(c);
    if fluid.adsInxSft == 2
        ads = max(ads, adsmax);
    end

end
