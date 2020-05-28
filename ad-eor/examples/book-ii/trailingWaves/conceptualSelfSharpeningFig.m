%% Conceptual illustration of the balance between self sharpening and smearing
% The parameter Sx defines how nonlinear the flux function will be compared
% with the f(u)=u.
set(gcf,'Position',[400 480 556 260]);
u = ones(1,20); u(11:20)=0;
x = 0.5:1:19.5;
s = linspace(0,1,51);
clf
[I,J,Sx,Sy] = deal(2:19, 8:13, .5, 10);
for i=1:14
    u(I) = u(I) + .1*(u(I-1)+u(I+1)-2*u(I));
end
cols = lines(6);
plot([.5 10 10 19.5],Sy*[1 1 0 0],'--k','LineWidth',1);
hold on
plot(x,Sy*u,'-',x(J),Sy*u(J),'o','LineWidth',2,'MarkerSize',8,...
    'Color',cols(1,:),'MarkerFaceColor',[.7 .7 1]);
vx = 0.5-u; vy=0*vx;
quiver(x(J),Sy*u(J),vx(J),vy(J),.5,'LineWidth',2,'Color',cols(2,:));
vx = 2*u - 1;
quiver(x(J),Sy*u(J),Sx*vx(J),vy(J),.5*Sx,'LineWidth',2,'Color',cols(5,:));
hold off
axis([-.5 20.5 -2 12]);
h = get(gca,'Children');
legend(h(1:2),'Self-sharpening','Smearing')
axis off
set(gca,'FontSize',12)

axes('Position',[.16 .19 .18 .5]);
f = (s.^2-s);
df = 2*s-1;
plot(s,Sx*f+s,s,0*s+s,'LineWidth',1);
legend('h(u)', 'u','Location','northwest','EdgeColor','none','Color','none')
set(gca,'FontSize',12,'XTick',[],'YTick',[]);
