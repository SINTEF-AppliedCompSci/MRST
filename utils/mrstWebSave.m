function downloadFcn = mrstWebSave()
%Get Call-Back for Downloading Online Resources Specified by URLs
%
% SYNOPSIS:
%   downloadFcn = mrstWebSave()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   downloadFcn - Call-Back function (function handle) for downloading
%                 online resources specified through a URL.  In recent
%                 versions of MATLAB this is just function `websave`. We
%                 wrap function `urlwrite` as a backwards compatbility
%                 measure in earlier versions of MATLAB.
%
% NOTE:
%   This function is mainly intended to support an encompassing download
%   manager such as function `githubDownload`.
%
%   The call-back function supports the following syntax
%
%      file = downloadFcn(file, url)
%      file = downloadFcn(file, url, 'p1', v1, ...)
%
%   in which the parameters are interpreted as follows
%
%      file - Name of local file into which contents of remote resource
%      will be saved.
%
%      url  - Uniform resource locator of remote resource (file contents).
%
%      'pn'/pv - List of key/value pairs that will be passed through to the
%      underlying GET request of the URL.  Numeric arguments converted to
%      strings using function `num2str`.
%
%   and the return value is
%
%      file - Unmodified input file name if successful, empty in case of
%      download failure.
%
% SEE ALSO:
%   `githubDownload`.

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

   if exist('websave', 'file') == 2
      downloadFcn = @websave;
   else
      downloadFcn = @web_save;
   end
end

%--------------------------------------------------------------------------

function file = web_save(file, url, varargin)
%Restricted fall-back option for WEBSAVE on older releases of MATLAB

   args = stringify_values(varargin{:});

   if ~isempty(args)
      args = [ { 'Get' }, { args } ];
   end

   [file, ok] = urlwrite(url, file, args{:});

   if ~ ok
      % Failed to download file from URL.
      file = [];
   end
end

%--------------------------------------------------------------------------

function args = stringify_values(varargin)
   vals = varargin(2 : 2 : end);

   i = cellfun(@(x) isnumeric(x) || islogical(x), vals);
   strng = cellfun(@num2str, vals(i), 'UniformOutput', false);

   [vals{i}] = strng{:};

   args = reshape([ varargin(1 : 2 : end) ; vals ], 1, []);
end
