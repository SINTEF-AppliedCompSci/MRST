function b = VEM2D_locb(X, hK, xK, yK, vol, PNstar)
        

%%  CELL DATA                                                            %%

n = size(X,1);

                                    %   Monomial products. mm(i) = m_i*m_j.
mm = @(X) 1/hK.*[ X(:,1)-xK,            X(:,2)-yK                  , ... 
                 (X(:,1)-xK).^2./hK,   (X(:,1)-xK).*(X(:,2)-yK)./hK , ...
                                        (X(:,2)-yK).^2./hK              ];

I = polygonInt(X,mm);
H = zeros(3,3);
H(1,:) = [vol, I(1:2)]; H(2,2:3) = I(3:4); H(3,3) = I(5);
H = H + tril(H',-1);
C = H*PNstar(1:3,1:n);
b = 1
