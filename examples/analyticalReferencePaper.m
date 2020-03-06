function [u, f, mu] = analyticalReferencePaper(d, kappa)
    
% d = spatial dimension
%
% Reference solution from paper
%  title={Finite volume methods for elasticity with weak symmetry},
%  author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
%  journal={International Journal for Numerical Methods in Engineering},
%  volume={112},
%  number={8},
%  pages={939--962},
%  year={2017},
%  publisher={Wiley Online Library}
%
%  The displacement u which is returned by this function has not been
%  multiplied by the indicator function, as it appears in the paper. 


    % displacement
    u = cell(d, 1);
    % force
    f = cell(d, 1);

    switch d 
      case 2
        
        mu = @(x, y) (1 - Xindfunction2D(x, y)) + kappa*Xindfunction2D(x, y);
        u{1} = @(x, y) 1./mu(x, y).*((x - 0.5).^2.*(y - 0.5).^2);
        u{2} = @(x, y) 1./mu(x, y).*(-2/3*(x - 0.5).*(y - 0.5).^3);
        
        f{1} = @(x, y) 0.25*(2*x - 1).^2 + 0.25*(2*y - 1).^2;
        f{2} = @(x, y) -2.0*x.*y + 1.0*x + 1.0*y - 0.5;
        
      case 3

        mu = @(x, y) (1 - Xindfunction3D(x, y)) + kappa*Xindfunction3D(x, y);
        u{1} = @(x,y,z) 1./mu(x, y, z).*(((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2);
        u{2} = @(x,y,z) 1./mu(x, y, z).*(((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2);
        u{3} = @(x,y,z) 1./mu(x, y, z).*(-2/3*((x - 0.5).*(y - 0.5).^2 + (x - 0.5).^2.*(y - 0.5)).*(z - 0.5).^3);

        f{1} = @(x,y,z) (0.0625*(2*x - 1).^2.*(2*z - 1).^2 + 1/4*(2*x - 1).*(y ...
                                                          - 0.5).*(2*z - 1).^2 ...
                         + 0.125*(2*x - 1).*(2*y - 1).^2.*(2*z - 1) - 1/6*(2*x ...
                                                          - 1).*(2*y - 1).*(3.0*z ...
                                                          - 1.5).*(x + y - 1) ...
                         + 0.125*(2*y - 1).^2.*(2*z - 1).^2);

        f{2} = @(x,y,z)((x - 0.5).*(2*y - 1).*(2*z - 1).^2/4 + 0.125*(2*x - ...
                                                          1).^2.*(2*y - 1).*(2*z ...
                                                          - 1) + 0.125*(2*x - ...
                                                          1).^2.*(2*z - 1).^2 ...
                        - (2*x - 1).*(2*y - 1).*(3.0*z - 1.5).*(x + y - 1)/6 ...
                        + 0.0625*(2*y - 1).^2.*(2*z - 1).^2);

        f{3} =  @(x,y,z) (0.0625*(2*x - 1).^2.*(2*z - 1).^2 - (2*x - 1).*(2*y ...
                                                          - 1).*(6.0*z - ...
                                                          3.0).*(x + y - 1)/6 ...
                          + (2*x - 1).*(2*z - 1).^2.*(-1.5*x - 3.0*y + 2.25)/12 ...
                          + 0.0625*(2*y - 1).^2.*(2*z - 1).^2 + (2*y - 1).*(2*z ...
                                                          - 1).^2.*(-3.0*x - ...
                                                          1.5*y + 2.25)/12); ...
                
      otherwise
        error('d not recognized')

    end
    
    
    
end

function u = Xindfunction2D(x, y)
    ix = (x > 0.5) & (y > 0.5);
    u = 0*x;
    u(ix) = 1;
end

function u = Xindfunction3D(x, y, z)
    ix = (x > 0.5) & (y > 0.5) & (z > 0.5);
    u = 0*x;
    u(ix) = 1;
end