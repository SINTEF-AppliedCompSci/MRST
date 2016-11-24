function u=cuw(u0, cfl, dx, T, flux, dflux, boundary, lim)
if nargin==7,
    lim = @(r) max(0,min(1.3*r,min(.5+.5*r,1.3)));
end
limiter = @(a,b) lim(b./(a+1e-6)).*a;
dt = cfl*dx/max(abs(dflux(u0)));

[u,U] = deal(u0); f=0*u;
t = 0.0;
n = numel(u0);
i = 3:n-2;
j = 1:n-1;
while (t<T)
  if (t+dt > T), dt = (T-t); end
  t=t+dt;

  u = boundary(u,2);
  du = u(2:end)-u(1:end-1);
  phi = limiter(du(1:end-1),du(2:end));
  phi = phi([end 1:end 1]);     % only correct for periodic b.c.
  ur = u + .5*phi;  fr = flux(ur); dfr = dflux(ur);
  ul = u - .5*phi;  fl = flux(ul); dfl = dflux(ul);
  ap = max(max(dfr,dfl),0);
  am = min(min(dfr,dfl),0);
  mdf = max(ap,-am); mdf=max(mdf);
  f(j) = (ap(j).*fr(j) - am(j+1).*fl(j+1) ...
      + ap(j).*am(j+1).*(ul(j+1) - ur(j)))./ (ap(j)-am(j+1)+1e-6);
  U(i) = u(i) - dt/dx*(f(i)-f(i-1));
  
  U = boundary(U,2);
  du = U(2:end)-U(1:end-1);
  phi = limiter(du(1:end-1),du(2:end));
  phi = phi([end 1:end 1]);     % only correct for periodic b.c.
  ur = U + .5*phi;  fr = flux(ur); dfr = dflux(ur);
  ul = U - .5*phi;  fl = flux(ul); dfl = dflux(ul);
  ap = max(max(dfr,dfl),0);
  am = min(min(dfr,dfl),0);
  mdf = max(max(ap,-am),mdf); mdf=max(mdf);
  f(j) = (ap(j).*fr(j) - am(j+1).*fl(j+1) ...
      + ap(j).*am(j+1).*(ul(j+1) - ur(j)))./ (ap(j)-am(j+1)+1e-6);
  u(i) = .5*u(i) +.5*( U(i)-dt/dx*(f(i)-f(i-1)) );
  dt = cfl*dx/mdf;
end