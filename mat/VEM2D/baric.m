function [E,BB,w]=baric(G);
% given G matrix defining anti-clockwise vertexes of local element ...
% ... gives area E and baricenter BB of element ...
% ... and w weight vector useful for quadrature at nodes (exact for linears)
% ... \int_E f dx \sim \sum_{j=1}^m w(j) f(node_j)   
% -------------------------------------------------------------------
m=size(G,1); % number of vertezes
xx=sum(G)/m;
BB=zeros(1,2);
Et=zeros(m,1);
BBt=zeros(m,2);
w=zeros(m,1);  % vector of weights (useful for load integration formula)
ww=0;
for j=1:m   % subdivide in sub-triangles and compute area and baricenter of each
    M=[G(j,1:2);G(mod(j,m)+1,1:2);xx];
    Et(j)=abs(det( [M,ones(3,1)] ))/2;  % area triangle 
    BBt(j,1:2)= sum(M)/3;  % baricenter triangle
    w(j) = w(j) + Et(j)/3;
    w(mod(j,m)+1) = w(mod(j,m)+1) + Et(j)/3;
    ww = ww + Et(j)/3;
end
E=sum(Et);
BB=(Et'*BBt)/E;
w = w + ones(m,1)*(ww/m);


% ----- CHECK option for triangular case ------------
%figure(1)
%hold on
%plot(0,0)
%plot(1.2,1.2)
%plot(G(1,1),G(1,2),'x')
%plot(G(2,1),G(2,2),'x')
%plot(G(3,1),G(3,2),'x')
%plot(BB(1),BB(2),'rx')
% ---------------------------------------------------