function area=heronsformula(vertex1,vertex2,vertex3,tol)
% Heron's numerically stable formula for calculating the area of a triangle
% defined by vertices in 3D space. Vertex is a 1x3 vector containing the
% x,y,z coordinates of the vertices.

a=norm(vertex1-vertex2);
b=norm(vertex1-vertex3);
c=norm(vertex2-vertex3);
rootarg=-det([0,a^2,b^2,1;a^2,0,c^2,1;b^2,c^2,0,1;1,1,1,0]);
if rootarg<tol
    rootarg=0;
end
area=0.25*sqrt(rootarg);

end
