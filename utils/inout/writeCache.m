function writeCache(arg, varargin)
% Write callers workspace to matfile in directory ./.cache/
%
% SYNOPSIS:
%   writeCache(fn)
%   writeCache(fn, {'var1', 'var2', ...})
%   writeCache({arg1, ...}, {'var1', 'var2', ...})
%
% PARAMETERS:
%   All variables in callers workspace.
%
% RETURNS:
%   none
%
% SEE ALSO:
%   `readCache`.

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

opt = struct('verbose', mrstVerbose);
 if nargin == 2
   var = varargin{1};
 else
   var = {};
 end

 if iscell(arg)
   if ~isdir(fullfile(pwd, '.cache'))
     mkdir('.cache');
   end
   id = md5sum(arg);
   fn = ['.cache', filesep, id];

   %dispif(true, 'Saving to %s\n', fn);
   dispif(opt.verbose, '--------------------------------------------------\n');
   dispif(opt.verbose, 'Writing');
   if opt.verbose, fprintf(' ''%s'',', var{:});end
   dispif(opt.verbose, '\b to file\n%s.mat\n', fn);
   dispif(opt.verbose, '--------------------------------------------------\n');

   evalin('caller', ['save ', fn, sprintf(' %s', var{:})]);

 else

 [path, name, type]=fileparts(arg);


 % old hash exist and is equal to current hash
 fp = fopen(['.cache', filesep, name, type], 'r');
 h  = hash(arg);

 if ~isdir('.cache')
   mkdir('.cache');

 elseif fp ~= -1
   oldhash = fgets(fp);

   %if hash has not changed
   if strcmp(oldhash, h)
     dispif(true, 'Not writing cache file.');
     return

   else
     dispif(true, 'Removing   .cache/%s.mat\n', oldhash);
   %if hash has changed
     system(['rm .cache', filesep, oldhash,'.mat']);
     system(['rm .cache', filesep, name, type]);
   end
 end

 fp = fopen(['.cache', filesep, name, type], 'w');
 fprintf(fp, '%s', h);
 fclose(fp);

 % Save
 dispif(opt.verbose, '--------------------------------------------------\n');
 dispif(opt.verbose, 'Writing ');
 if opt.verbose, fprintf('''%s'', ', var{:});end
 dispif(opt.verbose, ['to .cache', filesep, '%s.mat\n'], h);
 dispif(opt.verbose, '--------------------------------------------------\n');
 evalin('caller', ['save .cache', filesep, h, sprintf(' %s',var{:})]);
end
