%% Introduction to Automatic Differentiation
% Automatic differentiation is, quoting from
% <https://en.wikipedia.org/wiki/Automatic_differentiation Wikipedia>: "_a
% set of techniques to numerically evaluate the derivative of a function
% specified by a computer program. AD exploits the fact that every computer
% program, no matter how complicated, executes a sequence of elementary
% arithmetic operations (addition, subtraction, multiplication, division,
% etc.) and elementary functions (exp, log, sin, cos, etc.). By applying
% the chain rule repeatedly to these operations, derivatives of arbitrary
% order can be computed automatically, accurately to working precision, and
% using at most a small constant factor more arithmetic operations than the
% original program._"
%
% In MRST, we have implemented automatic differentiation using operator
% overloading. In this tutorial, we give you a very brief introduction to
% our implementation of AD, which is specifically geared towards
% discretization of partial differential equations and computation with
% long vectors and sparse matrices. Unlike many other implementations of
% AD, we have chosen to keep individual Jacobians as a list of sparse
% matrices and not automatically expand into a large global matrix. This
% ensures efficiency and flexibility in manipulating matrix blocks
% corresponding to subequations of (complex) equation systems.
%
% Suggested reading:
% 
% # R. D. Neidinger. Introduction to automatic differentiation and MATLAB
% object-oriented programming. SIAM Rev., 52(3):545–563, 2010. DOI:
% 10.1137/080743627
% # K.-A. Lie. An introduction to reservoir simulation using MATLAB: User
% guide for the Matlab Reservoir Simulation Toolbox (MRST). SINTEF ICT,
% December 2015, http://www.sintef.no/Projectweb/MRST/Publications


%% Simple formula using two scalars
% Compute z = 3 exp(-xy) and its partial derivatives wrt x and y. First, we
% compute it manually and then using automatic differentiation (AD). To
% understand what is going on behind the curtains when using AD, you can
% set a breakpoint at line 13 and use the 'Step into' function to step
% through the different functions that are invoked when evaluating z
[x,y]=deal(1,2);
disp('Manually computed:')
disp([1 -y -x]*3*exp(-x*y))

disp('Automatic differentiation:')
[x,y] = initVariablesADI(1,2);
z = 3*exp(-x.*y)                                                           %#ok<*NOPTS>


%% Demonstrate the importance of vectorization
m = 14;
[n,t1,t2,t3,t4] = deal(zeros(m,1));
clc
disp('Table: timing for differentiation of z=xy');
disp('   nvec  analytic     vectorized   loop+mtimes  loop+times');
disp('   --------------------------------------------------------');
for i = 1:m
   n(i) = 2^(i-1);
   xv = rand(n(i),1); yv=rand(n(i),1);
   [x,y] = initVariablesADI(xv,yv);
   tic, z = xv.*yv; zx=yv; zy = xv; t1(i)=toc;                             %#ok<*NASGU>
   tic, z = x.*y;                   t2(i)=toc;
   if i<17
      tic, for k =1:n(i), z(k)=x(k)*y(k); end;  t3(i)=toc;
      tic, for k =1:n(i), z(k)=x(k).*y(k); end; t4(i)=toc;
   end
   fprintf('%7d: %6.5f sec, %6.5f sec, %6.5f sec, %6.5f sec\n', ...
      n(i), t1(i), t2(i), t3(i), t4(i))
end
disp('   --------------------------------------------------------');
disp('');
loglog(n,t1,'-*',n,t2,'-+',n,t3,'-o',n,t4,'-s');
legend('analytical','vectorized','for+mtimes','for+times','Location','NorthEast');


%% Use automatic differentation to assmble a linear system
x = initVariablesADI(zeros(3,1));
eq1 = [ 3,  2, -4]*x + 5;
eq2 = [ 1, -4,  2]*x + 1;
eq3 = [-2,- 2.  4]*x - 6;
eq = cat(eq1,eq2,eq3);
disp(eq.jac{1})
disp(eq.val)
u  = -eq.jac{1}\eq.val; disp(u)

%% Use AD to solve the nonlinear Rosenbrock problem
% Problem: minimize f(x,y) = (a-x)^2 + b(y-x^2)^2 = 0
% The solution is (a,a^2). A necessary condition for minimum is that
% grad(f(x,y))=0, which gives two residual equations
%    2(a-x) - 4bx(y-x^2) = 0
%              2b(y-x^2) = 0
% Here, we set a=1 and b=100
[a,b,tol] = deal(1,100,1e-6);
[x0,incr] = deal([-.5;4]);
while norm(incr)>tol
    x = initVariablesADI(x0);
    eq = cat( 2*(a-x(1)) - 4*b.*x(1).*(x(2)-x(1).^2), ...
        2*b.*(x(2)-x(1).^2));
    incr = - eq.jac{1}\eq.val;
    x0 = x0 + incr;
end

%% Slightly more elaborate version with plotting
[a,b] = deal(1,100);
fx = @(x) 2*(a-x(1)) - 4*b.*x(1).*(x(2)-x(1).^2);
fy = @(x) 2*b.*(x(2)-x(1).^2);
[tol,res,nit] = deal(1e-6, inf,0);
x = [-.5 4];
while res>tol && nit<100
    nit = nit+1;
    xd = initVariablesADI(x(end,:)');
    eq = cat(fx(xd),fy(xd));
    incr = - eq.jac{1}\eq.val;
    res = norm(incr);
    x = [x; x(end,:)+incr']; %#ok<AGROW>
end

% Plot solution: use modified colormap to clearly distinguish the region
% with lowest colors. Likewise, show the value of f per iteration on a
% nonlinear scale, i.e., f(x,y)^0.1
subplot(4,1,1:3);
f = @(x,y) (a-x).^2 + b.*(y-x.^2).^2;
[xg, yg] = meshgrid(-2:0.2:2, -1.5:0.2:4);
surf(xg,yg,reshape(f(xg(:), yg(:)),size(xg))); shading interp
hold on; 
plot3(x(:,1),x(:,2),f(x(:,1),x(:,2)),'-or', 'LineWidth', 2, ...
    'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r'); 
hold off
view(-130,45); axis tight; colormap([.3 .3 1; jet(99)]);
subplot(4,1,4)
bar(f(x(:,1),x(:,2)).^.1,'b'); axis tight; set(gca,'YTick',[]);



%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
