function c = getVEColors(ind)
% Return colors used for plotting VE models
%
% SYNOPSIS:
%   c = getVEColors(ind)
%
% PARAMETERS:
%   ind - a scalar or vector consisting of one or more of the numbers one to
%         five. These numbers signify the following parts of a VE model:
%            1 - trapped CO2
%            2 - CO2 plume (movable)
%            3 - residual CO2
%            4 - CO2 dissolved in brine
%            5 - brine
%
% RETURNS:
%   c  - a numel(ind)x3 vector specifying one RGB color per entry in 'ind'
%
% EXAMPLE:
%   image(repmat(1:5,2,1)); colormap(getVEColors(1:5))
%
% SEE ALSO:
%   getInventoryColors

red = [1 0 0]; blue = [0 0 1]; white = [1 1 1]; black = [0 0 0];
mix = @(a,b,c) (a*c + b*(100-c))/100;

c = zeros(numel(ind),3);
for i=1:numel(ind)
   switch ind(i)
      case 1
         c(i,:) = mix(mix(red,white,40),black,80);
      case 2
         c(i,:) = mix(red,white,40);
      case 3
         c(i,:) = mix(mix(red,blue,50),white,20);
      case 4
         c(i,:) = mix(blue,white,20);
      case 5
         c(i,:) = mix(blue,white,30);
   end
end 