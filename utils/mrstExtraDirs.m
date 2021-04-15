function p = mrstExtraDirs(varargin)
%Get List of Directories Added to MATLAB's Search PATH by MRST
%
% SYNOPSIS:
%   p = mrstExtraDirs()
%   p = mrstExtraDirs(search)
%
% PARAMETERS:
%   search - String or cell-array of strings containing partial paths that
%            will be matched against the full list of MRST's extra
%            directories.  Only those directories that match at least one
%            of the search strings will be returned.  Nested cell-array
%            structure is not preserved.
%
%            OPTIONAL: If empty or not present, all of MRST's extra
%            directories will be returned.
%
%            String matching is aware of platform conventions.  We use case
%            insensitive matching on Microsoft Windows, and case sensistive
%            matching on Unix-like systems.
%
% RETURNS:
%   p - Cell array of strings with each element being a path to a directory
%       added by activating MRST or through the MRST module system.
%
% EXAMPLES:
%   % 1) Show all MRST directories
%   strvcat(mrstExtraDirs())
%
%   % 2) Show all MRST directories that match the string (partial path)
%   %       FULLFILE('dataset_manager', 'datasets')
%   strvcat(mrstExtraDirs(fullfile('dataset_manager', 'datasets')))
%
%   % 3) Show all MRST directories that match either of the strings
%   %       'wells_and_bc', 'testgrids'
%   %
%   % a) As separate string arguments.
%   strvcat(mrstExtraDirs('wells_and_bc', 'testgrids'))
%
%   % b) As a single cell-array of strings.
%   strvcat(mrstExtraDirs({'wells_and_bc', 'testgrids'}))
%
%   % c) As a single cell-array with nested cell-arrays of strings.
%   strvcat(mrstExtraDirs({'wells_and_bc', {'testgrids'}}))
%
%   % d) As a mix of strings and multiply nested cell-arrays of strings.
%   strvcat(mrstExtraDirs('wells_and_bc', {{{'testgrids'}}}))
%
% SEE ALSO:
%   ROOTDIR, mrstPath, mrstModule, addpath.

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

   if ~exist('ROOTDIR', 'file')
      error('MRST:NotActive', ...
           ['Function ''%s'' supported only in the ', ...
            'context of an active MRST session'], mfilename());
   end

   p = regexp(path, escape(pathsep), 'split');
   p = exclude_matlabroot(p);
   p = exclude_editordir(p);
   p = translate_rootdir(p);
   p = translate_tempdir(p);
   p = translate_homedir(p);

   if nargin > 0
      p = dir_subset(p, nonempty(varargin));
   end
end

%--------------------------------------------------------------------------

function s = escape(s)
   s = regexptranslate('escape', s);
end

%--------------------------------------------------------------------------

function s = dirstring(s)
   s = escape(canonicalise_dirname(s));
end

%--------------------------------------------------------------------------

function dname = canonicalise_dirname(dname)
   w     = what(fileparts(fullfile(dname, '.')));
   dname = w.path;
end

%--------------------------------------------------------------------------

function p = exclude_matlabroot(p)
   p = p(~ strncmp(matlabroot, p, numel(matlabroot)));
end

%--------------------------------------------------------------------------

function p = exclude_editordir(p)
   tdir = escape(fullfile(tempdir(), 'x'));

   m = regexp(p, ['^', tdir(1 : end - 1), '[Ee]ditor_.*'], 'match');
   p = p(cellfun('isempty', m));
end

%--------------------------------------------------------------------------

function p = translate_rootdir(p)
   for k = 0 : 2
      up  = repmat({'..'}, [1, k]);
      src = fullfile( ROOTDIR , up{:});
      dst = fullfile('ROOTDIR', up{:});

      p = regexprep(p, dirstring(src), escape(dst));
   end
end

%--------------------------------------------------------------------------

function p = translate_tempdir(p)
   p = regexprep(p, dirstring(tempdir()), 'tempdir');
end

%--------------------------------------------------------------------------

function p = translate_homedir(p)
   if ispc
      home = 'USERPROFILE';
   else
      home = 'HOME';
   end

   p = regexprep(p, dirstring(getenv(home)), 'HOME');
end

%--------------------------------------------------------------------------

function arr = nonempty(arr)
   arr = arr(~cellfun(@isempty, arr));
end

%--------------------------------------------------------------------------

function p = dir_subset(p, search)
   m = matcher();

   if isempty(search)
      i = true(size(p));
   else
      while any(cellfun(@iscell, search))
         % Flatten nested structures.
         search = [ search{:} ];
      end

      i = false(size(p));

      for s = reshape(search, 1, [])
         if isempty(s{1})
            % An empty string/array matches everything.  There's no need to
            % search any further because we've found our fixpoint.
            i(:) = true;
            break
         end

         i = i | ~cellfun(@isempty, m(p, s{1}));
      end
   end

   p = p(i);
end

%--------------------------------------------------------------------------

function m = matcher()
   if ispc()
      m = @(txt, patt) strfind(lower(txt), lower(patt));
   else
      m = @strfind;
   end
end
