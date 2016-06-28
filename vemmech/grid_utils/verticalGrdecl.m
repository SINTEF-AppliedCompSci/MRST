function grdecl = verticalGrdecl(grdecl, varargin)
    opt = struct('method', 'all')
    opt = merge_options(opt, varargin{:});
    top = [1, 2];
    bot = [4, 5];
    [xyz, zcorn] = grdeclXYZ(grdecl);
    switch opt.method
      case 'sides'
        xyz(bot, 1, :) = xyz(top, 1, :); % x-sides
        xyz(bot, end, :) = xyz(top, end, :);
        xyz(bot, :, 1) = xyz(top, :, 1); % y-sides
        xyz(bot, :, end) = xyz(top, :, end)
      case 'all'
        xyz(top, :, :) = xyz(bot, :, :); 
      otherwise
        error('no such method');     
    end
    grdecl.COORD = xyz(:);
end