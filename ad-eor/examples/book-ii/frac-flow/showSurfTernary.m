%% Illustrate Riemann solution for ternary displacement
% This script shows the graphical construction of the Riemann solution for
% a setup in which a surfactant solution is continuously injected into a
% reservoir that has been preflushed by water. 
%
% NB! This is not a general Riemann solver and will thus not work for some
% parameter choices. In these cases, you should first start by finding the
% point sw where the backward secant of the fw curve from s=sw to s=sor
% touches the fs curve. If we call intersection point ss, you then find the
% upper envelope of fs between s=ss and s=sors. 
col = lines(6);

%% Compute and plot the fractional flow curve for the right and left states
[swr, sor, ads] = deal(0.3,0.7,0);
[nw,no,M] = deal(2,3,2);
[swrs,sors,Ms] = deal(.05,.975,1);
s  = linspace(0,1,501);
iors = find(s<=sors,1,'last');
iwrs = find(s<=swrs,1,'last');
ior  = find(s<=sor,1,'last');
iwr  = find(s<=swr,1,'last');

sm = (min(max(s,swr),sor)-swr)/(sor-swr);
fw = sm.^nw./(sm.^nw + M*(1-sm).^no);
sm = (min(max(s,swrs),sors)-swrs)/(sors - swrs);
fs = sm.^nw./(sm.^nw + Ms*M*1*(1-sm).^no);

subplot(1,2,1),cla, hold on
set(gca,'FontSize',12)
plot(s,fw,'-','LineWidth',2,'Color',col(1,:));
plot(s,fs,'-','LineWidth',2,'Color',col(2,:));
plot(s([iwr ior]), fw([iwr ior]),'+','LineWidth',2,'Color',col(1,:));
plot(s([iwrs iors]), fs([iwrs iors]),'+','LineWidth',2,'Color',col(2,:));

%% Find envelope for surfactant curve
i3 = find(fs(2:iors)./(s(2:iors)+ads)-diff(fs(2:iors+1))./diff(s(2:iors+1))>0,1);
plot([-ads s(i3:iors)],fs([1 i3:iors]),'--','LineWidth',2,'Color',col(3,:));
plot(s(i3),fs(i3),'o','Color',col(3,:),'MarkerFaceColor',col(3,:));

%% Find intersection of this envelope and the pure water curve
i2 = find(min(fs(i3)/(s(i3)+ads).*(s(2:end)+ads),1)<fw(2:end),1);
plot(s(i2),fw(i2),'o','Color',col(5,:),'MarkerFaceColor',col(5,:));

%% Find envelope for pure water curve
flag = fw(i2:ior)./(s(i2:ior)-sor)>diff(fw(i2:ior+1))./diff(s(i2:ior+1));
i1 = min(find(flag)+i2,ior);
if isempty(i1), i1=ior; end
plot(s(i1(1)),fw(i1(1)),'s','Color',col(4,:),'MarkerFaceColor',col(5,:));
plot(s([i2 i1:ior]),fw([i2 i1:ior]),'--','LineWidth',2,'Color',col(5,:))

%% Plot self-similar solution
subplot(1,2,2)
hold all
xi = fliplr(diff(fs(i3:iors+1))./diff(s(i3:iors+1)));
ss = fliplr(s(i3:iors));
df = fliplr(diff(fw([i2 i1:ior]))./diff(s([i2 i1:ior])));
ds = s([i2 i1:ior]);
xi = [xi([1:end end]), df([1 1:end])];
ss = [ss, ds([1 1:end])];
xi = [xi([1:end end]), max(xi)+.5];
ss = [ss([1:end end]), sor];
plot(xi(ss<=sors),ss(ss<=sors),'-','LineWidth',2);
axis tight, set(gca,'Ylim',[0 1],'FontSize',12);
text(.05,sors,'S_{or}^s','HorizontalAlignment','left',...
    'VerticalAlignment','middle','FontSize',12);
text(.5*(xi(end-1)+xi(end)),sor,'1-S_{or}','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',12);
set(gca,'YLim',[min(ss)-.1 1]);

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
