function [CXXFLAGS, LINK, LIBS] = mrstDefaultMexFlags(defines)
%Define Common Compiler and Linker Flags/Libraries for MRST's MEX Functions
%
% SYNOPSIS:
%   CXXFLAGS               = mrstDefaultMexFlags()
%   CXXFLAGS               = mrstDefaultMexFlags(defines)
%   [CXXFLAGS, LINK]       = mrstDefaultMexFlags(...)
%   [CXXFLAGS, LINK, LIBS] = mrstDefaultMexFlags(...)
%
% PARAMETERS:
%   defines - Character vector, or cell-array of character vectors, of
%             user-specified preprocessor symbols that will be included in
%             the list of compiler flags and passed on to the compiler as
%             '#define'-d symbols.  OPTIONAL.  No additional preprocessor
%             symbols defined if this parameter is not specified.
%
% RETURNS:
%   CXXFLAGS - Compiler flags, including flags to enable OpenMP if
%              available in the primary toolchain.  Assumed to be forwarded
%              on to a C++ compiler.  Cell-array of character vectors.
%
%   LINK     - Additional flags to pass on to the toolchain linker.
%              Typically includes the path to MATLAB's binary/DLL/DSO
%              directory for the current architecture (COMPUTER('arch')).
%              Cell-array of character vectors.
%
%   LIBS     - Additional run-time libraries that are typically needed in
%              MRST's MEX functions.  Includes references to MATLAB's
%              built-in LAPACK and BLAS implementations (Intel's MKL) and
%              typically also an OpenMP runtime support library.
%              Cell-array of character vectors.
%
% NOTE:
%   This function assumes that the MRST MEX function in turn is implemented
%   in the C++ language.
%
% EXAMPLE:
%   % Define set of compiler/linker/library flags for a MEX function whose
%   % implementation can optionally use OpenMP and the AMGCL package if
%   % activated at compile time through preprocessor symbols.
%
%   [CXXFLAGS, LINK, LIBS] = ...
%       mrstDefaultMexFlags({'USE_OPENMP', 'USE_AMGCL'})
%
% SEE ALSO:
%   `buildmex`, `mex.getCompilerConfigurations`, `computer`.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

   if nargin == 0,     defines = {}; end
   if ischar(defines), defines = {defines}; end

   [LINK, mwlib, iomp5] = link_libraries();

   if is_visual_cpp() || is_visual_cpp_intel()
      [CXXFLAGS, iomp5] = compile_flags_msvc(defines, iomp5);

   elseif is_xcode_clang()
      [CXXFLAGS, iomp5] = compile_flags_clang(defines);

   elseif is_gnu_gcc()
      [CXXFLAGS, LINK, iomp5] = compile_flags_gcc(defines, LINK, iomp5);

   else
      report_error_unsupported_archictecture();
   end

   LIBS = [ iomp5, { mwlib('lapack'), mwlib('blas') } ];
end

%--------------------------------------------------------------------------

function [LINK, mwlib, iomp5] = link_libraries()
   arch = computer('arch');

   if ispc()
      mwlib = @(lib) ...
         fullfile(matlabroot, 'extern', 'lib', arch, ...
                  'microsoft', ['libmw', lib, '.lib']);

      LINK  = { ['-L', fullfile(matlabroot, 'bin', arch) ] };
      iomp5 = { 'libiomp5md.lib' };

   else
      mwlib = @(lib) ['-lmw', lib];
      LINK  = { ['-L', fullfile(matlabroot, 'sys', 'os', arch)] };
      iomp5 = { '-liomp5' };
   end
end

%--------------------------------------------------------------------------

function [CXXFLAGS, iomp5] = compile_flags_msvc(defines, iomp5)
   % Note explicit /EHsc to enable C++ exception handling
   CXXFLAGS = { ['COMPFLAGS=/EHsc /MD ', formatDefs('/', defines), ...
                 ' /openmp /wd4715 /fp:fast /O2 /bigobj'] };
   iomp5    = { ['LINKFLAGS=$LINKFLAGS ', ...
                 '/nodefaultlib:vcomp ', iomp5{1} ]};
end

%--------------------------------------------------------------------------

function [CXXFLAGS, iomp5] = compile_flags_clang(defines)
   dispif(mrstVerbose(), 'Clang detected. Will not use OpenMP...\n');

   CXXFLAGS = ...
      { ['CXXFLAGS=$CXXFLAGS ', formatDefs('-', defines), ...
         ' -fPIC -O3 -std=c++11 -ffast-math -march=native'] };

   iomp5 = {};
end

%--------------------------------------------------------------------------

function [CXXFLAGS, LINK, iomp5] = compile_flags_gcc(defines, LINK, iomp5)
   compile_native = '-march=native';
   if ispc()
      compile_native = '';
   end

   CXXFLAGS = ...
      { ['CXXFLAGS=$CXXFLAGS -D_GNU_SOURCE ', ...
         formatDefs('-', defines), ' -fPIC -O3 -std=c++11 ', ...
         '-ffast-math ', compile_native, ' -fopenmp'] };

   if ispc() || gcc_need_builtin_openmp()
      % Could be MinGW or similar on Windows or a GCC version that's too
      % new to integrate nicely with MATLAB's bundled copy of Intel's
      % OpenMP runtime libraries.  Use compiler's OpenMP runtime library.
      LINK  = [ LINK, 'LDFLAGS="$LDFLAGS -fopenmp"' ];
      iomp5 = {};
   end
end

%--------------------------------------------------------------------------

function report_error_unsupported_archictecture()
   error('Architecture:Unsupported', ...
        ['Computer Architecture ''%s'' (compiler ''%s'') is not ', ...
         'Supported for %s'], computer(), compilername(), mfilename());
end

%--------------------------------------------------------------------------

function s = formatDefs(prefix, defines)
   if isempty(defines)
      s = '';
   else
      s = sprintf([' ', prefix, 'D%s'], defines{:});
   end
end

%--------------------------------------------------------------------------

function tf = is_visual_cpp()
   tf = ispc() && ~isempty(regexpi(compiler_short_name(), '^MSVC'));
end

%--------------------------------------------------------------------------

function tf = is_visual_cpp_intel()
   tf = ispc() && ...
      ~isempty(regexpi(compiler_short_name(), 'INTELCPP\d+MSVCPP\d+'));
end

%--------------------------------------------------------------------------

function tf = is_gnu_gcc()
   tf = ~isempty(regexpi(compiler_short_name(), 'g\+\+'));
end

%--------------------------------------------------------------------------

function tf = is_xcode_clang()
   tf = ~isempty(regexpi(compiler_short_name(), 'Clang\+\+'));
end

%--------------------------------------------------------------------------

function tf = gcc_need_builtin_openmp()
   gccver = gcc_major_version();

   tf = ~all(isfinite(gccver)) || (gccver(1) > 7);
end

%--------------------------------------------------------------------------

function ver = gcc_major_version()
   cfg = compiler_config();

   if ~isempty(strtrim(cfg.Version))
      % Compiler version provided by configuration.  String of the form
      %
      %    7
      %
      % or
      %
      %    9.2.1
      %
      % Parse out the major/main component of this string.
      ver = parse_gcc_verstring(cfg.Version);
   else
      % Empty version string.  Fall back to increasingly desperate attempts
      % to parse output from the compiler itself.  Return 'NaN' if unable
      % to determine version.
      ver = parse_gcc_version_output(cfg);
   end
end

%--------------------------------------------------------------------------

function ver = parse_gcc_version_output(cfg)
%First fall-back strategy for GCC version if ISEMPTY(config.Version)
%
% Attempts to use
%
%   g++ -dumpversion
%
% as the source of version information.

   [stat, verstr] = system([cfg.Location, ' -dumpversion']);

   if stat == 0
      % g++ -dumpversion successful.  Output (verstr) is of the easy form
      %
      %   5.5.0
      %
      % or
      %
      %   7
      %
      % Convert this to numeric format.
      ver = parse_gcc_verstring(verstr);
   else
      % g++ -dumpversion failed.  Retry using "g++ --version" (more brittle
      % and expensive).
      ver = parse_gcc_fullversion_output(cfg);
   end
end

%--------------------------------------------------------------------------

function ver = parse_gcc_fullversion_output(cfg)
%Second fall-back strategy for GCC version if ISEMPTY(config.Version)
%
% Attempts to use
%
%   g++ --version
%
% as the source of version information.  Leads to a somewhat brittle and
% expensive computation.

   [stat, verstr] = system([cfg.Location, ' --version']);

   if stat ~= 0
      % g++ --version failed.  Can't infer version number.
      ver = NaN;
   else
      % g++ --version successful.  Output (verstr) is expected to be a text
      % block of the form
      %
      %   g++ (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
      %   Copyright (C) 2017 Free Software Foundation, Inc.
      %   This is free software; see the source for copying conditions.  There is NO
      %   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
      %
      % Attempt to extract and parse the number from the 'g++ ...' line.
      name   = regexptranslate('escape', cfg.ShortName);
      patt   = '\d+\.\d\+\.\d+';
      verstr = regexp(verstr, '\n', 'split');

      isver  = ~cellfun(@isempty, regexp(verstr, [name, '.*', patt]));
      verstr = regexp(verstr{isver}, [name, '.*(', patt, ')'], 'tokens');

      ver = parse_gcc_verstring(verstr{1}{1});
   end
end

%--------------------------------------------------------------------------

function ver = parse_gcc_verstring(verstr)
   ver = sscanf(strtrim(verstr), '%d.%d.%d');

   if numel(ver) < 3
      % E.g. verstr = '7'.  Zero-fill to three components.
      ver(3) = 0;
   end
end

%--------------------------------------------------------------------------

function cname = compilername()
   cfg   = compiler_config();
   cname = cfg.Name;
end

%--------------------------------------------------------------------------

function cname = compiler_short_name()
   cfg   = compiler_config();
   cname = cfg.ShortName;
end

%--------------------------------------------------------------------------

function cfg = compiler_config()
   cfg = mex.getCompilerConfigurations('c++');
end
