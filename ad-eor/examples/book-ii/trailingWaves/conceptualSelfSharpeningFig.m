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

%% Copyright Notice
%
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
