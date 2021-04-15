function tf = mrstIsLiveEditorDir(dirname)
%Detect if Script is Run From Live Editor or in Cell Mode
%
% SYNOPSIS:
%   tf = mrstIsLiveEditorDir(dirname)
%
% PARAMETERS:
%   dirname - Directory name.  Character vector or string if availble in
%             the current version of MATLAB.  Usually derived from a
%             statement of the form::
%
%                 dirname = fileparts(mfilename('fullpath'))
%
% RETURNS:
%   tf - Logical flag identifying if the `dirname` is likely to be a Live
%        Editor directory or if the `dirname` is the containing directory
%        of a script that is run interactively in Cell Mode.  TRUE in case
%        of live editor directory or empty `dirname`, FALSE otherwise.
%
% SEE ALSO:
%   `fileparts`, `mfilename`.

%{
Copyright 2020-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   tf = isempty(dirname) || is_live_editor_dir(dirname);
end

%--------------------------------------------------------------------------

function tf = is_live_editor_dir(dirname)
   edlivedir = regexptranslate('escape', fullfile(tempdir(), 'Editor'));

   tf = ~isempty(regexp(dirname, ['^', edlivedir], 'once'));
end
