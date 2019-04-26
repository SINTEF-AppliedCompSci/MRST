function [val, bz, mat] = readMatrixMarket(filename,is_matrix)
if(is_matrix)
    fid = fopen(filename);
    lin= fgetl(fid);
    lin= fgetl(fid);
    lin = split(lin)
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
    lin = split(lin)
    bz(1)= str2num(lin{4});
    bz(2)= str2num(lin{5});
    assert(str2num(lin{5})==1);
    sss = fscanf(fid,'%d',2)';
    rhs = fscanf(fid,'%f', sss(1));
    fclose(fid)
    val = rhs
    mat= val;
end
