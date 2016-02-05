
names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));

res = cell(numel(names),2);
for i = 1:numel(names)
    
    fprintf('--- Formation: %s -------\n', names{i});
    [~, ~, coarsening] = cell_size_coarsening_factors( names{i}, 2500, 16 );
    
    % store results:
    res{i,1} = names{i};
    res{i,2} = coarsening;

end

% save res to .mat file



% %%% Using coarsening_limit of 20, cell_size = 3000 meters:

% --- Formation: Tubaenfm -------
% Loading module libgeometry
% Max coarsening level possible is   20.
% Coarsening of  6 gives    605 cells a   3000 meter cell size.
% --- Formation: Stofm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives    605 cells a   3000 meter cell size.
% --- Formation: Bjarmelandfm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives   3284 cells a   3000 meter cell size.
% --- Formation: Arefm -------
% Max coarsening level possible is   20.
% Coarsening of 15 gives   2061 cells a   3000 meter cell size.
% --- Formation: Ilefm -------
% Max coarsening level possible is   20.
% Coarsening of 15 gives   2264 cells a   3000 meter cell size.
% --- Formation: Garnfm -------
% Max coarsening level possible is   20.
% Coarsening of 15 gives   1855 cells a   3000 meter cell size.
% --- Formation: Tiljefm -------
% Max coarsening level possible is   20.
% Coarsening of 15 gives   2166 cells a   3000 meter cell size.
% --- Formation: Brentgrp -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives   2108 cells a   3000 meter cell size.
% --- Formation: Brynefm -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives   4761 cells a   3000 meter cell size.
% --- Formation: Sleipnerfm -------
% Max coarsening level possible is   10.
% Coarsening of  3 gives    414 cells a   3000 meter cell size.
% --- Formation: Sognefjordfm -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives    908 cells a   3000 meter cell size.
% --- Formation: Fensfjordfm -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives    887 cells a   3000 meter cell size.
% --- Formation: Krossfjordfm -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives    919 cells a   3000 meter cell size.
% --- Formation: Huginfmeast -------
% Max coarsening level possible is   11.
% Coarsening of  3 gives    185 cells a   3000 meter cell size.
% --- Formation: Huginfmwest -------
% Max coarsening level possible is   10.
% Coarsening of  3 gives    330 cells a   3000 meter cell size.
% --- Formation: Sandnesfm -------
% Max coarsening level possible is   20.
% Coarsening of  3 gives   4628 cells a   3000 meter cell size.
% --- Formation: Ulafm -------
% Max coarsening level possible is   10.
% Coarsening of  3 gives    360 cells a   3000 meter cell size.
% --- Formation: Gassumfm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives    692 cells a   3000 meter cell size.
% --- Formation: Johansenfm -------
% Max coarsening level possible is   20.
% Coarsening of 15 gives    257 cells a   3000 meter cell size.
% --- Formation: Pliocenesand -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives    284 cells a   3000 meter cell size.
% --- Formation: Skadefm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives   1173 cells a   3000 meter cell size.
% --- Formation: Statfjordfm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives   3117 cells a   3000 meter cell size.
% --- Formation: Utsirafm -------
% Max coarsening level possible is   20.
% Coarsening of  6 gives   2280 cells a   3000 meter cell size.

%%% for 2500 meters
% --- Formation: Tubaenfm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives    908 cells a   2500 meter cell size.
% --- Formation: Stofm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives    908 cells a   2500 meter cell size.
% --- Formation: Bjarmelandfm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives   4845 cells a   2500 meter cell size.
% --- Formation: Arefm -------
% Max coarsening level possible is   15.
% Coarsening of 12 gives   3298 cells a   2400 meter cell size.
% --- Formation: Ilefm -------
% Max coarsening level possible is   15.
% Coarsening of 12 gives   3618 cells a   2400 meter cell size.
% --- Formation: Garnfm -------
% Max coarsening level possible is   15.
% Coarsening of 12 gives   2961 cells a   2400 meter cell size.
% --- Formation: Tiljefm -------
% Max coarsening level possible is   15.
% Coarsening of 12 gives   3456 cells a   2400 meter cell size.
% --- Formation: Brentgrp -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives   5019 cells a   2000 meter cell size.
% --- Formation: Brynefm -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives  11197 cells a   2000 meter cell size.
% --- Formation: Sleipnerfm -------
% Max coarsening level possible is   10.
% Coarsening of  2 gives   1131 cells a   2000 meter cell size.
% --- Formation: Sognefjordfm -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives   2194 cells a   2000 meter cell size.
% --- Formation: Fensfjordfm -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives   2159 cells a   2000 meter cell size.
% --- Formation: Krossfjordfm -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives   2235 cells a   2000 meter cell size.
% --- Formation: Huginfmeast -------
% Max coarsening level possible is   11.
% Coarsening of  2 gives    491 cells a   2000 meter cell size.
% --- Formation: Huginfmwest -------
% Max coarsening level possible is   10.
% Coarsening of  2 gives   1140 cells a   2000 meter cell size.
% --- Formation: Sandnesfm -------
% Max coarsening level possible is   15.
% Coarsening of  2 gives  10840 cells a   2000 meter cell size.
% --- Formation: Ulafm -------
% Max coarsening level possible is   10.
% Coarsening of  2 gives    967 cells a   2000 meter cell size.
% --- Formation: Gassumfm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives   1068 cells a   2500 meter cell size.
% --- Formation: Johansenfm -------
% Max coarsening level possible is   15.
% Coarsening of 12 gives    430 cells a   2400 meter cell size.
% --- Formation: Pliocenesand -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives    437 cells a   2500 meter cell size.
% --- Formation: Skadefm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives   1772 cells a   2500 meter cell size.
% --- Formation: Statfjordfm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives   4558 cells a   2500 meter cell size.
% --- Formation: Utsirafm -------
% Max coarsening level possible is   15.
% Coarsening of  5 gives   3422 cells a   2500 meter cell size.