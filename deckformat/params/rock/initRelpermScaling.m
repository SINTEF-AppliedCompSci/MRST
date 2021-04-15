function krscale = initRelpermScaling(deck, nc)
% set up the relperm scaling parameters from deck ('SOGCR', 'SGU', 'SGCR',
% 'SWCR', 'SOWCR', 'KRW', ....)
% Note: non-active cells are also assigned a value.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    drain = struct();
    imb = struct();
    [drain.w, drain.ow, drain.og, drain.g, ok_d] = getThreePhaseScaling(deck, ...
                                                      nc, '');
    [imb.w, imb.ow, imb.og, imb.g, ok_i] = getThreePhaseScaling(deck, nc, 'I');
    
    krscale.drainage = drain;
    krscale.imbibition = imb;
end

function [w, ow, og, g, present] = getThreePhaseScaling(deck, nc, prefix)
   [w, okw] = getRelPermScaling(deck, nc, prefix, 'W');
   [ow, okow] = getRelPermScaling(deck, nc, prefix, 'OW');
   [og, okog] = getRelPermScaling(deck, nc, prefix, 'OG');
   [g, okg] = getRelPermScaling(deck, nc, prefix, 'G');
   
   present = okw || okow || okog || okg;
end

function [pts, present] = getRelPermScaling(deck, nc, prefix, phase)

% Connate, critical, first s for max kr, max kr
   pts = repmat([NaN, NaN, NaN, NaN], nc, 1);

   connate = [prefix, 'S', phase, 'L'];
   crit = [prefix, 'S', phase, 'CR'];
   maxs = [prefix, 'S', phase, 'U'];
   if strcmpi(phase, 'ow') || strcmpi(phase, 'og')
      % OW or OG relperm has same maximum value
      maxv = [prefix, 'KRO'];
   else
      maxv = [prefix, 'KR', phase];
   end

   flds = {connate, crit, maxs, maxv};
   present = false;
   for i = 1:numel(flds)
       f = flds{i};
       if isfield(deck.PROPS, f)
           pts(:, i) = deck.PROPS.(f);
           present = true;
       end
   end
end
