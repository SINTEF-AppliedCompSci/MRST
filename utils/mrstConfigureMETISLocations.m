function varargout = mrstConfigureMETISLocations(varargin)
%Configure Installed Location of METIS Library and Associate Header Files
%
% SYNOPSIS:
%         mrstConfigureMETISLocations('pn1', 'pv1', ...)
%   cfg = mrstConfigureMETISLocations()
%   cfg = mrstConfigureMETISLocations('pn1', 'pv1', ...)
%
% PARAMETERS:
%   include  - Configure include directory.  Directory is expected to
%              contain the <metis.h> header directly.  Option name can be
%              abbreviated to 'inc' or 'i' (case insensitively)-mostly for
%              the case of interactive configuration in the Command Window.
%
%   library  - Configure library directory.  Option name can be abbreviated
%              to 'lib' or 'l' (ell, case insensitively)-mostly for the
%              case of interactive configuration in the Command Window.
%
%   binary   - Configure dynamic link directory-needed for DLLs on Windows.
%              Option name can be abbreviated to 'bin' or 'dll' (case
%              insensitively)-mostly for the case of interactive
%              configuration in the Command Window.
%
%   libmetis - Configure METIS library name.  Typically just 'metis', but
%              can be overridden if the installed version has some kind of
%              name decoration (e.g., to distinguish debug/profile/release).
%
% NOTE:
%   If at all possible, we recommend linking to a static version of the
%   METIS library since this reduces the risk of conflicting with any
%   internal METIS DLL (SO) used by MATLAB.
%
% RETURNS:
%   cfg - METIS configuration structure.  Contains the following fields
%           * bindir   - METIS dynamic link library
%           * incdir   - METIS include directory (MEX -I)
%           * libdir   - METIS library directory (MEX -L)
%           * libmetis - METIS library name (MEX -l)
%
% SEE ALSO:
%   mrstDefaultMexFlags, buildmex, mex.

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

   persistent CFG

   if isempty(CFG)
      CFG = empty_configuration();
   end

   assert (mod(nargin, 2) == 0, 'Inputs must be key/value pairs');

   if ~ (isempty(varargin) || inputs_are_strings(varargin{:}))
      error('ConfigArgs:Unsupported', ...
            'Inputs must be sequence of character vectors');
   end

   if nargin > 0
      munlock
      CFG = assign_configuration(CFG, varargin{:});
      mlock
   end

   if nargout > 0
      varargout{1} = CFG;
   end
end

%--------------------------------------------------------------------------

function cfg = empty_configuration()
   cfg = struct('incdir', '', 'libdir', '', 'bindir', '', ...
                'libmetis', 'metis');
end

%--------------------------------------------------------------------------

function tf = inputs_are_strings(varargin)
   tf = iscellstr(varargin) || ...
         (exist('isstring', 'builtin') && ...
          (isstring(varargin) || is_mix_char_string(varargin{:})));
end

%--------------------------------------------------------------------------

function cfg = assign_configuration(cfg, varargin)
   dirtype = cellstr(varargin(1 : 2 : end));
   dirs    = cellstr(varargin(2 : 2 : end));

   for e = 1 : numel(dirtype)
      switch lower(dirtype{e})
         case 'vcpkg'
            cfg = vcpkg(cfg, dirs{e});

         case {'include', 'inc', 'i'}
            cfg = include(cfg, dirs{e});

         case {'library', 'lib', 'l'}
            cfg = library(cfg, dirs{e});

         case {'bin', 'dll', 'binary'}
            cfg = dynlink(cfg, dirs{e});

         otherwise
            warning('Dirtype:Unknown', ...
                    'Directory type ''%s'' is unknown.  Ignored', ...
                    dirtype{e});
      end
   end
end

%--------------------------------------------------------------------------

function cfg = vcpkg(cfg, pkgdir)
   root = fullfile(pkgdir, 'installed', vcpkg_triplet());

   cfg = include(cfg, fullfile(root, 'include'));
   cfg = library(cfg, fullfile(root, 'lib'));

   if ispc
      % Support .DLL case
      cfg = dynlink(cfg, fullfile(root, 'bin'));
   end

   cfg = metis(cfg, 'metis');
end

%--------------------------------------------------------------------------

function cfg = include(cfg, incdir)
   cfg.incdir = incdir;
end

%--------------------------------------------------------------------------

function cfg = library(cfg, libdir)
   cfg.libdir = libdir;
end

%--------------------------------------------------------------------------

function cfg = dynlink(cfg, bindir)
   cfg.bindir = bindir;
end

%--------------------------------------------------------------------------

function cfg = metis(cfg, libmetis)
   cfg.libmetis = libmetis;
end

%--------------------------------------------------------------------------

function tf = is_mix_char_string(varargin)
   tf = all(cellfun(@(s) ischar(s) || isstring(s), varargin));
end

%--------------------------------------------------------------------------

function triplet = vcpkg_triplet()
   switch computer('arch')
      case 'win64'  , triplet = 'x64-windows';
      case 'glnxa64', triplet = 'x64-linux';
      otherwise
         error('VCPKG:UnsupportedTriplet', ...
              ['Don''t know how to map Computer Architecture ''%s'' ', ...
               'to a VCPKG triplet\n', ...
               'Please use individual directory specifications'], ...
               computer('arch'));
   end
end
