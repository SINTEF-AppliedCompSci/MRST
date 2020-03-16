function [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags(defines)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
    if nargin == 0
        defines = {};
    end
    if ischar(defines)
        defines = {defines};
    end
    a = computer('arch');

    if ispc()
        mwlib = @(lib) ...
            fullfile(matlabroot, 'extern', 'lib', a, ...
            'microsoft', ['libmw', lib, '.lib']);

        LINK  = { ['-L', fullfile(matlabroot, 'bin', a) ] };
        iomp5 = { 'libiomp5md.lib' };
    else
        mwlib = @(lib) ['-lmw', lib];
        LINK  = { ['-L', fullfile(matlabroot, 'sys', 'os', a)] };
        iomp5 = { '-liomp5' };
    end

    if is_visual_cpp()
        % Note explicit /EHsc to enable C++ exception handling
        CXXFLAGS  = { [sprintf('COMPFLAGS=/EHsc /MD %s', formatDefs('-', defines)), ...
            '/openmp /wd4715 /fp:fast /O2 /bigobj'] };
        iomp5     = { ['LINKFLAGS=$LINKFLAGS ', ...
            '/nodefaultlib:vcomp ', iomp5{1} ]};
        libstdcpp = {};

    elseif is_xcode_clang()
        dispif(mrstVerbose(), 'Clang detected. Will not use OpenMP...\n');
        CXXFLAGS = ...
            { [sprintf('CXXFLAGS=$CXXFLAGS %s ', formatDefs('-', defines)), ...
            '-fPIC -O3 -std=c++11 -ffast-math -march=native'] };

        libstdcpp = {};
        iomp5 = {};

    elseif is_gnu_gcc()
        if ispc()
            march = '';
        else
            march = '-march=native';
        end
        CXXFLAGS = ...
            { [sprintf('CXXFLAGS=$CXXFLAGS -D_GNU_SOURCE %s ', formatDefs('-', defines)), ...
               sprintf('-fPIC -O3 -std=c++11 -ffast-math %s -fopenmp', march)] };

        libstdcpp = {};
        if ispc()
            LINK = [LINK, 'LDFLAGS="$LDFLAGS -fopenmp"'];
            libstdcpp = {};
            iomp5 = {};
        end
    else

        error('Architecture:Unsupported', ...
            ['Computer Architecture ''%s'' (compiler ''%s'') is ', ...
            'not Supported for %s'], ...
            computer(), compilername(), mfilename());

    end

    LIBS = [ iomp5, { mwlib('lapack'), mwlib('blas') }, libstdcpp ];
    end

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

function tf = is_gnu_gcc()
    tf = ~isempty(regexpi(compiler_short_name(), 'g\+\+'));
end

%--------------------------------------------------------------------------

function tf = is_xcode_clang()
    tf = ~isempty(regexpi(compiler_short_name(), 'Clang\+\+'));
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
