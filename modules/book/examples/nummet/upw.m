function u=upw(u0,cfl,dx,T,flux,df,boundary)
u = u0; t = 0.0;
dt = cfl*dx/max(abs(df(u0)));
i = 2:numel(u0)-1;
while (t<T)
  if (t+dt > T)
     dt = (T-t);
  end
  t=t+dt;
  u = boundary(u);
  f = flux(u);
  u(i) = u(i) - dt/dx*(f(i)-f(i-1));
  dt = cfl*dx/max(abs(df(u)));
end
