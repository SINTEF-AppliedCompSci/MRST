function v = applyFunction(fn, array, varargin)
%Apply a function to entries of (cell) arrays, giving cell arrays as output
%
% SYNOPSIS:
%   f = applyFunction(@(x) fn(x), array)
%   f = applyFunction(@(x, y) fn(x, y), array1, array2)
%
% REQUIRED PARAMETERS:
%   fn    - function_handle that will be applied to each element of all
%          subsequent inputs.
%
%   array - Either a cell array or a regular double array. The outputs will
%           match the dimensions of this input, as a cell array.
%
% OPTIONAL PARAMETERS:
%   extra - Any number of additional arrays that match array in type and
%           dimension can be passed as additional arguments.
%
% RETURNS:
%   v     - Cell array where each entry contains the function result. Even
%           when the function outputs are scalar, this will be a cell
%           array.
% EXAMPLE:
%   x = applyFunction(@(x) sum(x, 2), {rand(10, 5), rand(3, 2)})
%
% NOTE:
%   All calls are redirected to arrayfun or cellfun with the
%   UniformOutput parameter disabled. This function is primarily a
%   convenience function since this type of call occurs a lot in MRST, and
%   also allows for easier Octave compatibility.
%
% SEE ALSO:
%   arrayfun, cellfun

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

    if iscell(array)
        f = @cellfun;
    else
        f = @arrayfun;
    end
    v = f(fn, array, varargin{:}, 'UniformOutput', false);
end
