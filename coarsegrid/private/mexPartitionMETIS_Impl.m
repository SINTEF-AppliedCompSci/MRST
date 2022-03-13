function varargout = mexPartitionMETIS_Impl(varargin)
%MATLAB Gateway Routine to METIS Library Functions
%
% SYNOPSIS:
%   p = mexPartitionMETIS_Impl(G, numBlocks)
%   p = mexPartitionMETIS_Impl(G, numBlocks, opt)
%   p = mexPartitionMETIS_Impl(G, numBlocks, opt, w)
%
% DESCRIPTION:
%   Compiles and invokes MEX function of the same name.  Requires a C++11
%   compiler.
%
% PARAMETERS:
%   G         - Grid structure.
%
%   numBlocks - Number of blocks into which to partition the grid 'G'.
%
%   opt       - METIS options structure.
%
%   w         - Connection weighting array.  Size equal to number of
%               half-faces (size(G.cells.faces, 1), number of faces
%               (G.faces.num), or number of internal faces
%               (sum(all(G.faces.neighbors > 0, 2)).  Often just the
%               one-sided transmissibilities calculated by function
%               `computeTrans` or the total face transmissibilities
%               calculated by function `getFaceTransmissibility`.
%
% RETURNS:
%   p - Partition vector.  Contains no empty or multiply connected blocks.
%
% SEE ALSO:
%   `mexPartitionMETIS`.

%{
Copyright 2020-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   [CXXFLAGS, LINK, LIBS] = mrstDefaultMexFlags();

   cfg = mrstConfigureMETISLocations();
   src = {'mexPartitionMETIS_Impl.cpp'};
   opt = { '-O' };

   INC  = include_dir(cfg);
   LINK = library_dir(LINK, cfg);
   LIBS = [ LIBS, { ['-l', cfg.libmetis] } ];

   buildmex(opt{:}, INC{:}, CXXFLAGS{:}, src{:}, LINK{:}, LIBS{:});

   copy_dll_if_applicable(cfg);

   % Invoke actual MEX accelerator.
   [varargout{1 : nargout}] = mexPartitionMETIS_Impl(varargin{:});
end

%--------------------------------------------------------------------------

function INC = include_dir(cfg)
   INC = append_if_isdir({}, cfg.incdir, '-I');
end

%--------------------------------------------------------------------------

function LINK = library_dir(LINK, cfg)
   LINK = append_if_isdir(LINK, cfg.libdir, '-L');
end

%--------------------------------------------------------------------------

function copy_dll_if_applicable(cfg)
   metisdll = fullfile(cfg.bindir, dll_name(cfg));

   if ispc && ~isempty(cfg.bindir) && isdir(cfg.bindir) && ...
      exist(metisdll, 'file')
      [stat, msg] = copyfile(metisdll, fileparts(mfilename('fullpath')));

      if stat == 0
         warning('DLLCopy:Unsuccessful', ...
                ['Unable to copy METIS DLL to output directory: %s\n', ...
                 'MEX file may not run'], msg);
      end
   end
end

%--------------------------------------------------------------------------

function dlist = append_if_isdir(dlist, dname, optname)
   if ~isempty(dname) && isdir(dname)                           %#ok<ISDIR>
      dlist = [ dlist, { [optname, dname] } ];
   end
end

%--------------------------------------------------------------------------

function fname = dll_name(cfg)
   [fname, fname] = fileparts(cfg.libmetis);                    %#ok<ASGLU>

   fname = regexprep(fname, '^lib', '', 'ignorecase');
   fname = [fname, '.dll'];
end
