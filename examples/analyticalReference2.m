mrstModule add ad-core
clear all

N = 10;

xtbl.xind = (1 : N)';
xtbl = IndexTable(xtbl);

ytbl.yind = (1 : N)';
ytbl = IndexTable(ytbl);

ztbl.zind = (1 : N)';
ztbl = IndexTable(ztbl);

xyztbl = crossTable(xtbl, ytbl, {}, 'optpureproduct', true);
xyztbl = crossTable(xyztbl, ztbl, {}, 'optpureproduct', true);

x = linspace(0, 1, N)';
y = linspace(0, 1, N)';
z = linspace(0, 1, N)';

map = TensorMap();
map.fromTbl = xtbl;
map.toTbl = xyztbl;
map.mergefds =  {'xind'};
map = map.setup();
x = map.eval(x);

map = TensorMap();
map.fromTbl = ytbl;
map.toTbl = xyztbl;
map.mergefds =  {'yind'};
map = map.setup();
y = map.eval(y);

map = TensorMap();
map.fromTbl = ztbl;
map.toTbl = xyztbl;
map.mergefds =  {'zind'};
map = map.setup();
z = map.eval(z);


u1 = @(x,y,z) ((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2;
u2 = @(x,y,z) ((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2;
u3 = @(x,y,z) -2/3*((x - 0.5).*(y - 0.5).^2 + (x - 0.5).^2.*(y - 0.5)).*(z - 0.5).^3;


f1 = @(x,y,z) (0.0625*(2*x - 1)**2*(2*z - 1)**2 + (2*x - 1)*(1.0*y - 0.5)*(2*z ...
                                                  - 1)**2/4 + 0.125*(2*x - ...
                                                  1)*(2*y - 1)**2*(2*z - 1) - ...
               (2*x - 1)*(2*y - 1)*(3.0*z - 1.5)*(x + y - 1)/6 + 0.125*(2*y - ...
                                                  1)**2*(2*z - 1)**2);

f2 = @(x,y,z)((1.0*x - 0.5)*(2*y - 1)*(2*z - 1)**2/4 + 0.125*(2*x - 1)**2*(2*y ...
                                                  - 1)*(2*z - 1) + 0.125*(2*x ...
                                                  - 1)**2*(2*z - 1)**2 - (2*x ...
                                                  - 1)*(2*y - 1)*(3.0*z - ...
                                                  1.5)*(x + y - 1)/6 + 0.0625*(2*y ...
                                                  - 1)**2*(2*z - 1)**2);

f3 =  @(x,y,z) (0.0625*(2*x - 1)**2*(2*z - 1)**2 - (2*x - 1)*(2*y - 1)*(6.0*z ...
                                                  - 3.0)*(x + y - 1)/6 + (2*x ...
                                                  - 1)*(2*z - 1)**2*(-1.5*x - ...
                                                  3.0*y + 2.25)/12 + 0.0625*(2*y ...
                                                  - 1)**2*(2*z - 1)**2 + (2*y ...
                                                  - 1)*(2*z - 1)**2*(-3.0*x - ...
                                                  1.5*y + 2.25)/12);

u1d1 = @(x,y,z) (2*(x - 0.5).*(y - 0.5).^2.*(z - 0.5).^2);
u1d2 = @(x,y,z) (2*(x - 0.5).^2.*(y - 0.5).*(z - 0.5).^2);
u1d3 = @(x,y,z) (2*(x - 0.5).^2.*(y - 0.5).^2.*(z - 0.5));
u2d1 = @(x,y,z) (2*(x - 0.5).*(y - 0.5).^2.*(z - 0.5).^2);
u2d2 = @(x,y,z) (2*(x - 0.5).^2.*(y - 0.5).*(z - 0.5).^2);
u2d3 = @(x,y,z) (2*(x - 0.5).^2.*(y - 0.5).^2.*(z - 0.5));
u3d1 = @(x,y,z) (-2/3*(z - 0.5).^3.*((y - 0.5).^2 + 2*(x - 0.5).*(y - 0.5)));
u3d2 = @(x,y,z) (-2/3*(z - 0.5).^3.*(2*(x - 0.5).*(y - 0.5) + (x - 0.5).^2));
u3d3 = @(x,y,z) (-2*(z - 0.5).^2.*((x - 0.5).*(y - 0.5).^2 + (x - 0.5).^2.*(y - 0.5)));

s11 = @(x, y, z) (2*u1d1(x, y, z));
s22 = @(x, y, z) (2*u2d2(x, y, z));
s33 = @(x, y, z) (2*u3d3(x, y, z));
s12 = @(x, y, z) (u1d2(x, y, z) + u2d1(x, y, z));
s13 = @(x, y, z) (u1d3(x, y, z) + u3d1(x, y, z));
s23 = @(x, y, z) (u2d3(x, y, z) + u3d2(x, y, z));
s21 = s12;
s31 = s13;
s32 = s23;


dotest = false;
if dotest
    
    [xAD, yAD, zAD] = initVariablesADI(x, y, z);

    u1AD = u1(xAD, yAD, zAD);
    u2AD = u2(xAD, yAD, zAD);
    u3AD = u3(xAD, yAD, zAD);

    u1dxAD = diag(u1AD.jac{1});
    u1dyAD = diag(u1AD.jac{2});
    u1dzAD = diag(u1AD.jac{3});
    u2dxAD = diag(u2AD.jac{1});
    u2dyAD = diag(u2AD.jac{2});
    u2dzAD = diag(u2AD.jac{3});
    u3dxAD = diag(u3AD.jac{1});
    u3dyAD = diag(u3AD.jac{2});
    u3dzAD = diag(u3AD.jac{3});

    u1dxnum = u1d1(x, y, z);
    u1dynum = u1d2(x, y, z);
    u1dznum = u1d3(x, y, z);
    u2dxnum = u2d1(x, y, z);
    u2dynum = u2d2(x, y, z);
    u2dznum = u2d3(x, y, z);
    u3dxnum = u3d1(x, y, z);
    u3dynum = u3d2(x, y, z);
    u3dznum = u3d3(x, y, z);

    max(abs(u1dxnum - u1dxAD))
    max(abs(u1dynum - u1dyAD))
    max(abs(u1dznum - u1dzAD))
    max(abs(u2dxnum - u2dxAD))
    max(abs(u2dynum - u2dyAD))
    max(abs(u2dznum - u2dzAD))
    max(abs(u3dxnum - u3dxAD))
    max(abs(u3dynum - u3dyAD))
    max(abs(u3dznum - u3dzAD))
    max(abs(u3dznum - u3dzAD))

end