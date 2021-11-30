function varargout = md5sum(varargin)
% md5sum - Compute md5 check sum of all input arguments
%
% SYNOPSIS:
%   str = md5sum(args...)
%
% PARAMETERS:
%   vararign - `md5sum` can take an arbitrary number of arguments.
%
% RETURNS:
%   varargout - string of 32 characters with the hexadecimal md5 checksum
%               of all numeric and character arrays that are found as plain
%               arrays, in structs or cell arrays.
%
% EXAMPLE:
%   C{1}=struct('a',1,'b',2);C{3}=speye(4);
%   sum = md5sum(C)
%
% NOTE:
%   This utility is written in C.  It must be compiled with ::
%
%      mex md5sum.c
%
%   before use.  On older systems, the command is ::
%
%      mex -DOLDMATLAB md5sum.c

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

if ~isunix && usejava('jvm')
    % md5 sum fails to build on windows, fall back to java support.
    % If there is no java support, hope that the compiler miraculously
    % works...
    [varargout{1:nargout}] = md5sum_fallback(varargin{:});
    return
end

try
    % Build
    buildmex md5sum.c

    % run
    [varargout{1:nargout}] = md5sum(varargin{:});
catch ex
    if strcmp(ex.identifier, 'MATLAB:MEX:genericFailure')
        [varargout{1:nargout}] = md5sum_fallback(varargin{:});
    else
        rethrow(ex);
    end
end
