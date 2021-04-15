%%
% Test on command line:
x=linspace(0,1,100).';
y=rand(size(x));
xi=linspace(0,1,1e6).';

%%

N=10;

fprintf('----------------------------------------------\n');

t=tic;
for i=1:N
    interp1(x,y,xi);
end
t1=toc(t);

t=tic;
for i=1:N
    interp1q(x,y,xi);
end
t2=toc(t);

t=tic;
for i=1:N
    mex_nakeinterp1(x,y,xi);
end
t3=toc(t);


t1/t3
t2/t3


    %%

%mex -O -v nakeinterp1.c
mex -O -v CFLAGS="\$CFLAGS -std=c99" nakeinterp1.c

%%

xi=0.;
idx1 = nakeinterp1(x, y, xi)
idx2 = interp1(x, y, xi, 'linear', 'extrap')

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
