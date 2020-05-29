function deck = getDeckEGG(varargin)
% Get the parsed deck for the EGG model

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('realization', []);
    opt = merge_options(opt, varargin{:});
    
    mrstModule add deckformat ad-props
    % Read and process file.
    pth = getDatasetPath('egg');
    fn  = fullfile(pth, 'MRST', 'Egg_Model_ECL.DATA');

    deck = readEclipseDeck(fn);
    
    if ~isempty(opt.realization)
        fn  = fullfile(pth, 'Permeability_Realizations', ...
            ['PERM', num2str(opt.realization), '_ECL.INC']);
        nc = prod(deck.GRID.cartDims);
        fid = fopen(fn);
        kw = {'PERMX', 'PERMY', 'PERMZ'};
        deck.GRID = replaceKeywords(deck.GRID, fid, kw, nc);
    end
    deck = convertDeckUnits(deck);
end