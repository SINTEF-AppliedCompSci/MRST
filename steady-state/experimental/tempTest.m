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


