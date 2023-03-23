
%% Demonstrate various analytical fluid models
mrstModule add incomp

%% Simple fluid model
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
                       'rho', [1000, 1000]*kilogram/meter^3, ...
                       'n'  , [   3,   3]);
s=linspace(0,1,20)';
kr = fluid.relperm(s);
plot(s,kr(:,1),'-s',s,kr(:,2),'-o');


%% Corey fluid model
fluid = initCoreyFluid('mu' , [   1,  10]*centi*poise     , ...
                       'rho', [1014, 859]*kilogram/meter^3, ...
                       'n'  , [   3, 2.5]                 , ...
                       'sr' , [ 0.2, 0.15]                 , ...
                       'kwm', [   1, .85]);
s=linspace(0,1,50)';
[kr,dkr,ddkr] = fluid.relperm(s);
subplot(1,3,1), 
plot(s(s<.85),kr(s<.85,1),s(s>.2),kr(s>.2,2),'LineWidth',1); axis([0 1 -.02 1]);
subplot(1,3,2),
plot(s(s<.85),dkr(s<.85,1),s(s>.2),dkr(s>.2,4),'LineWidth',1); axis([0 1 -.1 4.5]);
subplot(1,3,3),
plot(s(s<.85),ddkr(s<.85,1),s(s>.2),ddkr(s>.2,2),'LineWidth',1); axis([0 1 -.2 9.5]);

%% Simple fluid with linear capillary term
fluid = initSimpleFluidPc('mu' , [   1,  10]*centi*poise     , ...
                          'rho', [1014, 859]*kilogram/meter^3, ...
                          'n'  , [   2,   2], ...
                          'pc_scale', 2*barsa);
s = linspace(0, 1, 101).'; kr = fluid.relperm(s);
subplot(1,2,1), plot(s, kr), legend('kr_1(S)', 'kr_2(S)')
x.s = [s 1-s]; pc = fluid.pc(x);
subplot(1,2,2), plot(s, pc); legend('P_c(S)');


%% Simple model with Leverett-J function
G = computeGeometry(cartGrid([40,5,1],[1000 500 1]));
state.s = G.cells.centroids(:,1)/1000;
p = G.cells.centroids(:,2)*.001;
rock = makeRock(G, p.^3.*(1e-5)^2./(0.81*72*(1-p).^2), p);
fluid = initSimpleFluidJfunc('mu' , [   1,  10]*centi*poise     , ...
            'rho', [1014, 859]*kilogram/meter^3, ...
            'n'  , [   2,   2], ...
            'surf_tension',10*barsa/sqrt(0.1/(100*milli*darcy)),...
            'rock',rock);
plot(state.s, fluid.pc(state)/barsa,'o');