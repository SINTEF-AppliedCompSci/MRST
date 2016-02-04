function u = periodic(u)
u(1) = u(end-1);
u(end) = u(2);