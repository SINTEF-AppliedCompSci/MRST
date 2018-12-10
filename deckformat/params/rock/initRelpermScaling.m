function krscale = initRelpermScaling(deck, nc)
% set up the relperm scaling parameters from deck ('SOGCR', 'SGU', 'SGCR',
% 'SWCR', 'SOWCR', 'KRW', ....)
% Note: non-active cells are also assigned a value.
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
   [g, okg] = getRelPermScaling(deck, nc, prefix, 'W');
   
   present = okw || okow || okog || okg;
end

function [pts, present] = getRelPermScaling(deck, nc, prefix, phase)

% Connate, critical, first s for max kr, max kr
   pts = repmat([NaN, NaN, NaN, NaN], nc, 1);

   connate = [prefix, 'S', phase, 'L'];
   crit = [prefix, 'S', phase, 'CR'];
   maxs = [prefix, 'S', phase, 'U'];
   if strcmpi(phase, 'wo') || strcmpi(phase, 'og')
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
