function AK = FEM2D_loc_sq(G,K, k)

addpath('../VEM2D/');

                            %   Lobatto shape function products
l0l0_x = @(X) ((1-X(:,1))/2).^2;
l1l0_x = @(X) (1+X(:,1))/2.*(1-X(:,1))/2;
l1l1_x = @(X) ((1+X(:,1))/2).^2;

l0l0_y = @(X) ((1-X(:,2))/2).^2;
l1l0_y = @(X) (1+X(:,2))/2.*(1-X(:,2))/2;
l1l1_y = @(X) ((1+X(:,2))/2).^2;




l0l0_x_int = polygonInt_v2(G, K, l0l0_x, k+1);
l1l0_x_int = polygonInt_v2(G, K, l1l0_x, k+1);
l1l1_x_int = polygonInt_v2(G, K, l1l1_x, k+1);
l0l0_y_int = polygonInt_v2(G, K, l0l0_y, k+1);
l1l0_y_int = polygonInt_v2(G, K, l1l0_y, k+1);
l1l1_y_int = polygonInt_v2(G, K, l1l1_y, k+1);

AK = zeros(4,4);
AK(1,1:4) = 1/4*[ l0l0_y_int + l0l0_x_int, ...
                     -l0l0_y_int + l1l0_x_int, ...
                     -l1l0_y_int - l1l0_x_int, ...
                      l1l0_y_int - l0l0_x_int];
AK(2,2:4) = 1/4*[ l0l0_y_int + l1l1_x_int, ...
                      l1l0_y_int - l1l1_x_int, ...
                     -l1l0_y_int - l1l0_x_int];
AK(3,3:4) = 1/4*[ l1l1_y_int + l1l1_x_int, ...
                     -l1l1_y_int + l1l0_x_int];
AK(4,4)   = 1/4*[ l1l1_y_int + l0l0_x_int];

AK = AK + tril(AK',-1);

end