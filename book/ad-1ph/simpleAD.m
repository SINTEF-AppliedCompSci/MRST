%% Simple formula using two scalars
[x,y] = initVariablesADI(1,2);
z = 3*exp(-x.*y)

%% Demonstrate the importance of vectorization
m = 15;
[n,t1,t2,t3,t4] = deal(zeros(m,1));
for i = 1:m
   n(i) = 2^(i-1);
   xv = rand(n(i),1); yv=rand(n(i),1);
   [x,y] = initVariablesADI(xv,yv);
   tic, z = xv.*yv; zx=yv; zy = xv; t1(i)=toc;
   tic, z = x.*y;                   t2(i)=toc;
   if i<17
      tic, for k =1:n(i), z(k)=x(k)*y(k); end;  t3(i)=toc;
      tic, for k =1:n(i), z(k)=x(k).*y(k); end; t4(i)=toc;
   end
   fprintf('%7d: %6.5f sec, %6.5f sec, %6.5f sec, %6.5f sec\n', ...
      n(i), t1(i), t2(i), t3(i), t4(i))
end
loglog(n,t1,'-*',n,t2,'-+',n,t3,'-o',n,t4,'-s');
legend('analytical','vectorized','for+mtimes','for+times',2);

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
