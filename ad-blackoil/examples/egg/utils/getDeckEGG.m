function deck = getDeckEGG(varargin)
%Undocumented Utility Function

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

    opt = struct('realization', 0);
    opt = merge_options(opt, varargin{:});
    
    mrstModule add deckformat ad-props
    % Read and process file.
    pth = getDatasetPath('egg');
    fn  = fullfile(pth, 'MRST', 'Egg_Model_ECL.DATA');

    deck = readEclipseDeck(fn);
    
    assert(opt.realization >= 0 && opt.realization <= 100, ...
        'Egg realization must be a number between 0 and 100.')
    if opt.realization > 0
        fn  = fullfile(pth, 'Permeability_Realizations', ...
            ['PERM', num2str(opt.realization), '_ECL.INC']);
        nc = prod(deck.GRID.cartDims);
        fid = fopen(fn);
        kw = {'PERMX', 'PERMY', 'PERMZ'};
        deck.GRID= replaceKeywords(deck.GRID, fid, kw, nc);
    end
    deck = convertDeckUnits(deck);
end
