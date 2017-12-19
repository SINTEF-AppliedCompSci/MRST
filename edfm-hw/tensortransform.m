function tensor_new=tensortransform(tensor_old,coords_old,coords_new)
% Transforms tensor_old from old coordinate system to new coordinate
% system. coords is a 3x3 matrix with each row being a UNIT vector pointing
% in 3 orthogonal directions

e1_o=coords_old(1,:); 
% e1_o=e1_o/norm(e1_o);
e2_o=coords_old(2,:); 
% e2_o=e2_o/norm(e2_o);
e3_o=coords_old(3,:); 
% e3_o=e3_o/norm(e3_o);

e1_n=coords_new(1,:); 
% e1_n=e1_n/norm(e1_n);
e2_n=coords_new(2,:); 
% e2_n=e2_n/norm(e2_n);
e3_n=coords_new(3,:); 
% e3_n=e3_n/norm(e3_n);

% http://www.continuummechanics.org/coordxforms.html
Q=[dot(e1_n,e1_o) dot(e1_n,e2_o) dot(e1_n,e3_o);
   dot(e2_n,e1_o) dot(e2_n,e2_o) dot(e2_n,e3_o);
   dot(e3_n,e1_o) dot(e3_n,e2_o) dot(e3_n,e3_o)];

tensor_new=Q*tensor_old*Q';


end