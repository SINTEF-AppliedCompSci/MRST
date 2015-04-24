function varargout = mex_compute_fs_flux(varargin)
%Derive fine-scale contact fluxes from basis-function projection
%
% SYNOPSIS:
%   flux = mex_compute_fs_flux(G, mex_cs, block_hflux)
%
% PARAMETERS:
%   G      - Fine-scale grid structure.
%
%   mex_cs - Coarse-system structure as defined by function
%            'mex_generate_coarsesystem'.  Assumed to correspond to the
%            fine-scale grid defined by 'G'.
%
%   block_hflux -
%            Coarse-scale half-fluxes.  Assumed to correspond to
%
%               faceFlux2cellFlux(CG, cellFlux2faceFlux(CG, v))
%
%            in which 'v' is the (first) return value from function
%            'mex_compute_press_flux' and 'CG' is the coarse grid
%            definition upon which the coarse-scale system is constructed.
%
% RETURNS:
%   flux - Fine-scale fluxes defined on the connections (faces) of the
%          fine-scale model 'G'.
%
% SEE ALSO:
%   mex_generate_coarsesystem, mex_compute_press_flux.

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

   INCLUDE = {['-I', fullfile(d, '..', 'include', 'opm', 'pressure')], ...
              ['-I', fullfile(d, '..', 'mrst_api')]};

   LINK = { ['-L', fullfile(d, '..', 'lib')] };

   OPTS = { '-g', '-largeArrayDims' };

   SRC = { 'mex_compute_fs_flux.c', 'mrst_msmfe_support.c', ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c') };

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_fs_flux(varargin{:});
end
