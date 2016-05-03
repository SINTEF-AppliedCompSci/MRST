%%
figure()
n=100;
[X,Y]=meshgrid(linspace(-2,2,n),linspace(-2,2,n));
i=complex(0,1);
Z=X+i*Y;
sq1 =@(x) sqrt(abs(x)).*exp(i*angle(x)/2)
sq2 =@(x) sqrt(abs(x)).*exp(i*(angle(x)+2*pi)/2)
%contourf(X,Y,log(-(sq1(Z-0.5)+sq2(Z+0.5))),10),colorbar
mesh(X,Y,real(log((sq1(Z-0.5)+sq2(Z+0.5)))))

%%
x=linspace(-10,10,1000);
fr1=@(Z) -(sq1(Z-0.5)+sq2(Z+0.5));
fr =@(Z) real(log(-(sq1(Z-0.5)+sq2(Z+0.5))))
figure(),plot(x,real(fr(x)),x,-log(2*sqrt(x)))

%%
x=linspace(-10,10,1000);
lw=0.3
ana=@(Z) -real(2*log(1/lw*(-(sq1(Z-lw)+sq2(Z+lw)))));
figure(),plot(x,ana(x),x,log(x))

%% 2d skize
figure()
n=100;
H=2;d=H*0.4;
[X,Y]=meshgrid(linspace(-4,4,n),linspace(-H,H,n));
i=complex(0,1);
Z=X+i*Y;
ww=@(z,d) (2*H/pi)*asinh( (exp((z*pi)/(2*H))+ i*sin(d*pi/(2*H))));
G=@(w,d) 2*H*log(2*(sinh((w*pi)/(2*H))-i*sin(d*pi/(2*H))))/pi;
%contourf(X,Y,log(-(sq1(Z-0.5)+sq2(Z+0.5))),10),colorbar
mesh(X,Y,real(G(Z,d)))
%%
d=linspace(-H,H,100);
plot(d,real(G(i*H,d)),d,0*d)

