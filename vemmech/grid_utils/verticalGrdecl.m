function grdecl = verticalGrdecl(grdecl, varargin)
%
%
% SYNOPSIS:
%   function grdecl = verticalGrdecl(grdecl, varargin)
%
% DESCRIPTION: Straightened up the pillars in a corner point grid, given in
% Eclipse, and make them vertical
%
% PARAMETERS:
%   grdecl   - Grid structure in Eclipse format
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   method - 'all' : all pillars are made vertical
%
%            'sides' : Only the pillare on the sides are made vertical
%
%
% RETURNS:
%   grdecl - Grid structure with vertical pillars in Eclipse format.
%
% EXAMPLE:
% SEE ALSO:
%

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