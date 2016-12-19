function deck = getDeckEGG(varargin)
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