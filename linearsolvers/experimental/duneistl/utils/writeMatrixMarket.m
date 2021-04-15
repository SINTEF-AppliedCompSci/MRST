function writeMatrixMarket(mat, bz, filename,is_matrix)
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
    [i,j,s] = find(mat);
    bmat=[i,j,s];
    fn = fopen(filename,'w');
    fprintf(fn,'%%%%MatrixMarket matrix coordinate real general\n');
    fprintf(fn,'%% ISTL_STRUCT blocked ');
    fprintf(fn,'%d ',bz(1));
    fprintf(fn,'%d',bz(2));
    fprintf(fn,'\n');
    msz = size(mat);
    fprintf(fn,'%d ',msz(1));
    fprintf(fn,'%d ',msz(2));
    fprintf(fn,'%d',size(bmat,1));
    fprintf(fn,'\n');
    [i,j,s] = find(mat);
    for k=1:size(bmat,1)
        fprintf(fn,'%d ',bmat(k,1));
        fprintf(fn,'%d ',bmat(k,2));
        fprintf(fn,'%g\n',bmat(k,3));
    end
%fprintf(fn,'\n');
   fclose(fn);
else
    assert(bz(2)==1);
    fn = fopen(filename,'w');
    fprintf(fn,'%%%%MatrixMarket matrix array real general\n');
    fprintf(fn,'%% ISTL_STRUCT blocked ');
    fprintf(fn,'%d ',bz(1));
    fprintf(fn,'%d',bz(2));
    fprintf(fn,'\n');
    msz = size(mat);
    fprintf(fn,'%d ',msz(1));
    fprintf(fn,'%d ',msz(2));
    %fprintf(fn,'%d',size(mat,1));
    fprintf(fn,'\n');
    fprintf(fn,'%g\n',mat);
    fclose(fn);
end
