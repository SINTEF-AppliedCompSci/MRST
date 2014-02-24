function c = getInventoryColors(ind)
% c = getInventoryColors(ind) returns a color map used for plotting CO2
% inventories using green to orange colors. The input argument 'ind' is a
% scalar or vector consisting of one or more of the numbers 1 to 7 that
% signify the following entities in the inventory:
%   1 - dissolved
%   2 - residual (traps)
%   3 - residual
%   4 - residual (plume)
%   5 - movable (traps)
%   6 - movable (plume)
%   7 - leaked
%
% To look at the colors
%   image(repmat(1:6,2,1)); colormap(getInventoryColors(2:7))
red    = [1 0 0];
green  = [0 1 0];
yellow = [1 1 0];
orange = [1 .65 0];
white  = [1 1 1]; 
black  = [0 0 0];
mix = @(a,b,c) (a*c + b*(100-c))/100;

c = zeros(numel(ind),3);
for i=1:numel(ind)
   switch ind(i)
      case 1
         c(i,:) = mix(green,black,50);
      case 2
         c(i,:) = mix(green,black,75);
      case 3
         c(i,:) = green;
      case 4
         c(i,:) = mix(green,white,50);
      case 5
         c(i,:) = mix(yellow,white,80);
      case 6
         c(i,:) = mix(orange,white,80);
      case 7
         c(i,:) = mix(orange,red,30);
   end
end 
