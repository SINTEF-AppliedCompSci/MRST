function u=nt(u0,cfl,dx,T,flux,dflux,boundary,lim)
if nargin==7,
    lim = @(r) max(0,min(1.3*r,min(.5+.5*r,1.3)));
end
limiter = @(a,b) lim(b./(a+1e-6)).*a;
dt = cfl*dx/max(abs(dflux(u0)));

u=u0; t = 0.0;
while (t<T)
  if (t+dt > T), dt = .5*(T-t); end

  u = boundary(u,2);
  du = u(2:end)-u(1:end-1);
  s = limiter(du(1:end-1),du(2:end));
  s = s([end 1:end 1]);
  f = flux(u);
  df = f(2:end)-f(1:end-1);
  sigma = limiter(df(1:end-1),df(2:end));
  sigma = sigma([end 1:end 1]);     % only correct for periodic b.c.
  v = u - .5*dt/dx*sigma;
  g = flux(v) + .125*dx/dt*s;
  u(2:end) = .5*(u(1:end-1)+u(2:end)) - dt/dx*(g(2:end)-g(1:end-1));

  u = boundary(u,2);
  du = u(2:end)-u(1:end-1);
  s = limiter(du(1:end-1),du(2:end));
  s = s([end 1:end 1]);
  f = flux(u);
  df = f(2:end)-f(1:end-1);
  sigma = limiter(df(1:end-1),df(2:end));
  sigma = sigma([end 1:end 1]);     % only correct for periodic b.c.
  v = u - .5*dt/dx*sigma;
  g = flux(v) + .125*dx/dt*s;
  u(1:end-1) = .5*(u(1:end-1)+u(2:end)) - dt/dx*(g(2:end)-g(1:end-1));

  t = t+2*dt;
  dt = cfl*dx/max(abs(dflux(u)));
end