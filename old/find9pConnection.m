function N=find9pConnection(Gtop)
% find neibours based on 9p connection
    ni1    = Gtop.faces.neighbors(:,1)+1;
    ni2    = Gtop.faces.neighbors(:,2)+1;
    matrix_con=sparse(double(ni1),double(ni2),1,Gtop.cells.num+1,Gtop.cells.num+1);
    matrix_con=matrix_con+matrix_con';
    matrix_con2=matrix_con*matrix_con;
    matrix_con=matrix_con+(matrix_con2==2);
    [ni1,ni2,s] = find(matrix_con);%#ok
    N=[ni1,ni2]-1;
    N=N(ni1>ni2,:);
end