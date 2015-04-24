function c = myCOColor(ind)
% c = myCOColor(ind) returns a color used for VE plotting where 'ind' is
% a scalar or vector consisting of one or more of the number 1 to 5. These
% numbers signify the following fluid type (and tikz coloring scheme)
%   1 - trapped   : red!40!white!80!black
%   2 - plume     : red!40
%   3 - residual  : red!50!blue!20
%   4 - dissolved : blue!20
%   5 - brine     : blue!30
%
% To look at the colors
%   image(repmat(1:5,2,1)); colormap(myCOColor(1:5))
red = [1 0 0]; blue = [0 0 1]; white = [1 1 1]; black = [0 0 0];
mix = @(a, b, c) (a * c + b * (100 - c)) / 100;

c = zeros(numel(ind), 3);
for i = 1:numel(ind)
   switch ind(i)
     case 1
       c(i, :) = mix(mix(red, white, 40), black, 80);
     case 2
       c(i, :) = mix(red, white, 40);
     case 3
       c(i, :) = mix(mix(red, blue, 50), white, 20);
     case 4
       c(i, :) = mix(blue, white, 20);
     case 5
       c(i, :) = mix(blue, white, 30);
   end
end