function success = readCache(arg, varargin)
%Read cached variables associated with given file into caller's workspace.
%
% SYNOPSIS:
%   success = readCache(filename)
%   success = readCache({variable1, ...})
%
% PARAMETERS:
%   file    - Name (string) of file with which to associate a cache.
%
% RETURNS:
%   success - Whether or not cache reading succeeded.  If successful,
%             the cached variables are loaded directly into the caller's
%             workspace.
%
% SEE ALSO:
%   `writeCache`.

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


opt  = struct('verbose', mrstVerbose);
opt  = merge_options(opt, varargin{:});
success = false;
if iscell(arg)

  id = md5sum(arg);
  fn = ['.cache', filesep, id, '.mat'];
  if exist(fn, 'file') == 2,
    dispif(opt.verbose, '--------------------------------------------------\n');
    dispif(opt.verbose, ['Reading ', fn, '\n']);

    stuff = load(fn);
    name  = fieldnames(stuff);

    verbose = opt.verbose && ~isempty(name);
    dispif(verbose, 'Retrieving');

    for i=1:numel(name),
       dispif(verbose, ' ''%s'',',name{i});
       assignin('caller', name{i}, stuff.(name{i}));
    end
    dispif(verbose, '\b.\n');
    dispif(opt.verbose, '--------------------------------------------------\n');

    success = true;

  end



else

  if ~ischar(arg),
    error('readCache takes a string or cell array ar argument');
  end


  [path, name, type] = fileparts(arg);
  fp = fopen(['.cache', filesep, name, type]);
  if fp > 0,
    % Retrieve old hash of arg from file .cache/fn
    oldhash = fgets(fp);
    fclose(fp);

    % Compare old hash to current hash
    if strcmp(hash(arg), oldhash)
      dispif(true, 'Reading .cache%s%s.mat\n', filesep, oldhash);
      evalin('caller', ['load .cache', filesep, oldhash]);
      success = true;
    end
  end
end
