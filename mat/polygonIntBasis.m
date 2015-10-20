function I = polygonIntBasis(X, Xmid, XB, vol, m)

Xdof = [X; Xmid; XB];

NK = size(Xdof,1);

tri = delaunay(Xdof);
nTri = size(tri,1);
                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
Xq = [0.0, 1.0; 0.0, 0.0];%; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
w = vol.*[1/36*ones(size(X,1),2); 1/9*ones(size(Xmid,1),2); [4/9, 4/9]];
                            %   For each triangle t, evaluate integral.
I = zeros(NK,2);
% for phii = 1:NK
%     triTmp = tri((sum(tri == phii,2) == 1),:)
% end
for t = 1:nTri
                            %   
    dofVec = tri(t,2:3);
                            %   Triangle points
    x1 = Xdof(tri(t,1),:); x2 = Xdof(tri(t,2),:); x3 = Xdof(tri(t,3),:);
                            %   map phi : tRef -> t
%     plot(Xdof([tri(t,:),tri(t,1)] ,1), Xdof([tri(t,:),tri(t,1)],2));
    psi = @(x) x*[x1(1)-x3(1)   x1(2)-x3(2);
                  x2(1)-x3(1)   x2(2)-x3(2)] + repmat(x3,size(x,1),1);
                            %   determinant of Dphi, Jacobian of phi.
    detDpsi = (x1(1)-x3(1))*(x2(2)-x3(2)) - (x1(2)-x3(2))*(x2(1)-x3(1));
                            %   Map quadrature points to t
%    Xdof(tri(t,:),:)
%     XIq = psi(Xq)
%     x2
%     x3
                            %   Evaluate integral
%     I(dofVec, :) = I(dofVec, :) + ...
%                          abs(detDpsi)*1/36.*m([x2;x3]);
end
I = w.*m(Xdof);