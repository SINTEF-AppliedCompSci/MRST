%% Evaluate the van Genuchten model
m=2/3;
mw = 1;
mn = 1;
s = linspace(0,1,80);
pc = (s.^(-1/m) - 1).^(1-m);
kw = s.^.5.*(1-(1-s.^(1/m)).^m).^2/mw;
kn = (1-s).^.5.*(1-s.^(1/m)).^(2*m)/mn;

%% Shown nonlinear functions entering formulation in phase pressures
% Originallly: P_c(S), k_rw(S) and k_rn(S)
% Transformed: S(p_c), k_rw(p_c), and k_rn(p_c)
args = {'Interpreter','latex','FontSize',12};
subplot(2,2,1);
plot(s, pc,'LineWidth',1);
h=legend('$P_c(S)$','Location','best'); set(h,args{:});
subplot(2,2,2);
plot(s,kw, s,kn,'LineWidth',1);
h=legend('$k_{rw}(S)$', '$k_{rn}(S)$','Location','best'); set(h,args{:});
subplot(2,2,3); 
plot(pc, s,'LineWidth',1); 
h=legend('$P_c^{-1}(p_c)$','Location','best'); set(h,args{:});
subplot(2,2,4); plot(pc,kw, pc, kn,'LineWidth',1);
h=legend('$k_{rw}(P_c^{-1}(p_c))$', '$k_{rn}(P_c^{-1}(p_c))$','Location','best'); set(h,args{:});

%% Show fractional flow functions
subplot(1,2,1),
plot(s, kw./(kw + kn), s, kw./ (kw + .2*kn), s, kw./(kw + 5*kn),'LineWidth',1);
h = legend('$\mu_n=\mu_w$', '$\mu_n=5\mu_w$', '$\mu_n=\mu_w/5$',2);
set(h,args{:});
subplot(1,2,2);
plot(s, kw.*kn./(kw + kn), s, kw.*kn./ (5*kw + kn), s, kw.*kn./(kw + 5*kn),'LineWidth',1);
h = legend('$\mu_n=\mu_w$', '$\mu_n=5\mu_w$', '$\mu_n=\mu_w/5$',2);
set(h,args{:});
