function [Csym,Casym] = splitStiffnessMatrix(C)
% Split elastic modulus, used in weakly symmetric version of mpsa.
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

Nc = numel(C);

Nd = sqrt(size(C{1},2));
assert(Nd == 2 || Nd == 3)

Csym = cell(Nc,1);
Casym = Csym;

for iter1 = 1 : Nc
    cl = C{iter1};
    s = diag(diag(cl));
    if Nd == 2
        s(1,4) = cl(1,4);
        s(4,1) = cl(4,1);
    else % 3D
        s(1,[5 9]) = cl(1,[5 9]);
        s(5,[1 9]) = cl(5,[1 9]);
        s(9,[1 5]) = cl(9,[1 5]);
    end
    a = cl - s;
    Csym{iter1} = s;
    Casym{iter1} = a;
end