function [a,b] = elasticModuloTransform(c,d,from,to)
%http://en.wikipedia.org/wiki/Shear_modulus
%http://en.wikipedia.org/wiki/Lam%C3%A9_parameters

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
if(strcmp(from,'E_nu') && strcmp(to,'lam_mu'))
    E=c;nu=d;
    lam=E.*nu./((1+nu).*(1-2*nu));
    mu=E./(2*(1+nu));
    a=lam;
    b=mu;
elseif  (strcmp(from,'lam_mu') && strcmp(to,'E_nu'))
    %lam G=mu
    lam=c;mu=d;
    nu=lam./(2*(lam+mu));
    E=mu.*(3*lam+2*mu)./(mu+lam);
    a=E;b=nu;
else
    error('Not implemented')
end

end

