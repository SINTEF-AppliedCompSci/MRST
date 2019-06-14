function writeMatrixMarket(mat, bz, filename,is_matrix)
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
