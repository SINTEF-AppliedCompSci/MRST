function varargout = mex_compute_coarse_BI(varargin)
%Compute coarse-scale inverse inner product from fine-scale contributions
%
% SYNOPSIS:
%   BI = mex_compute_coarse_BI(mex_cs, totmob)
%
% PARAMETERS:
%   mex_cs - Coarse-system structure as defined by function
%            'mex_generate_coarsesystem'.
%
%   totmob - Fine-scale total mobility field.  One (positive) scalar value
%            for each cell in the fine-scale grid model upon which the
%            coarse-system structure is defined.
%
% RETURNS:
%   BI - Coarse-scale inverse inner product values suitable for passing to
%        function 'mex_schur_comp_symm' for the purpose of computing
%        elemental contributions to the (coarse-scale) Schur-complement
%        system of linear equations
%
% SEE ALSO:
%   mex_generate_coarsesystem, mex_schur_comp_symm.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   d = fileparts(mfilename('fullpath'));

   CFLAGS = {'CFLAGS="\$CFLAGS', '-Wall', '-Wextra', '-ansi', ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align', ...
             '-Wpointer-arith', '-Wbad-function-cast', ...
             '-Wmissing-prototypes', '-Wstrict-prototypes', ...
             '-Wmissing-declarations', '-Winline', '-Wundef', ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow', ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   LDFLAGS = {'LDFLAGS="\$LDFLAGS', ...
              ['-Wl,-rpath=', fullfile(d, '..', 'lib'), '"']};

   INCLUDE = {['-I', fullfile(d, '..', 'include', 'opm', 'pressure')]};

   LINK = { ['-L', fullfile(d, '..', 'lib')] };

   OPTS = { '-O', '-largeArrayDims' };

   SRC = { 'mex_compute_coarse_BI.c', 'mrst_msmfe_support.c'};

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_coarse_BI(varargin{:});
end
