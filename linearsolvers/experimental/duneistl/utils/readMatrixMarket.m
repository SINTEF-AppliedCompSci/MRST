function [val, bz, mat] = readMatrixMarket(filename,is_matrix)
%Undocumented Utility Function

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

if(is_matrix)
    fid = fopen(filename);
    lin= fgetl(fid);
    lin= fgetl(fid);
    lin = split(lin);
    bz(1)= str2num(lin{4});
    bz(2)= str2num(lin{5});
    ss = fscanf(fid,'%d',3)';
    mat = fscanf(fid,'%d %d %f', [3,inf])';
    fclose(fid);
    istl_mat=sparse(mat(:,1),mat(:,2),mat(:,3),ss(1),ss(2));
    val = istl_mat;
else
    %%
    fid = fopen(filename);
    lin= fgetl(fid);
    lin= fgetl(fid);
    lin = split(lin);
    bz(1)= str2num(lin{4});
    bz(2)= str2num(lin{5});
    assert(str2num(lin{5})==1);
    sss = fscanf(fid,'%d',2)';
    rhs = fscanf(fid,'%f', sss(1));
    fclose(fid)
    val = rhs;
    mat= val;
end
