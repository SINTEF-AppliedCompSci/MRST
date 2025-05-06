function u=upwind(u0,dt,dx,T,flux,boundary)
u = u0;
t = 0.0; r = dt/dx;
i = 2:numel(u0)-1;
while (t<T)
  if (t+dt > T)
     dt = (T-t);
     r = dt/dx;
  end
  t=t+dt;
  u = boundary(u);
  f = flux(u);
  u(i) = u(i) - r*(f(i)-f(i-1));
end
