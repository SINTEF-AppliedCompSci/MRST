function W = setWellSign(W)
%Ensure that wells have a defined sign. Will attempt to guess based on controls.
%
% SYNOPSIS:
%   W = setWellSign(W)
%
% DESCRIPTION:
%   The AD based solvers assume a W.sign is defined. This routine attempts
%   to ensure that wells do have a sign. A positive sign is used to
%   indicate an injector and a negative sign for a producer. If the wells
%   have rate controls, they will be given signs based on the signs of the
%   rates. If they are pressure controlled wells, they will get sign 0.
%
% REQUIRED PARAMETERS:
%   W  - Well struct, from e.g. addWell.
%
% RETURNS:
%   W  - Well struct where numel(W(i).sign) is guaranteed to be 1.
%

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
    for i = 1:numel(W)
        if isempty(W(i).sign)
            warning('mrst:badWellSign', ...
                ['Well ', W(i).name, ' has empty sign, guessing type']);
            W(i).sign = determineSign(W(i));
        end
    end
end

function s = determineSign(w)
    if strcmpi(w.type, 'pressure') || strcmpi(w.type, 'bhp')
        s = 0;
    else
        s = 1 - 2*(w.val < 0);
    end
end
