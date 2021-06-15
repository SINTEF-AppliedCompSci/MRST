%% Illustrate Riemann solutions for continuous displacement
% This script shows the graphical construction of the Riemann solution for
% a setup in which either a diluted polymer solution or a surfactant is
% continuously injected into a reservoir filled with oil and water at
% connate saturation. In the polymer-flooding case, polymer and water are
% assumed to be fully mixing (i.e., we set w=1 in the Todd-Longstaff mixing
% model).
%
% The script has several adjustable parameters. To produce the polymer
% figure from the MRST book chapter with full mixing and no adsorption, set
% polymer to true and ads=0. For the polymer figure with adsorption, you
% should set ads=0.6 and rerun all sections of this script except for the
% last one.
polymer = false;

%% Compute and plot the fractional flow curve for the right and left states
if polymer
    [swr, sor, ads] = deal(0.2,0.8,0);
    % [swr, swor, ads] = deal(0.2,0.8,0.6);
    [nw,no,M] = deal(2,3,5);
    s  = linspace(0,1,501);
    sm = (min(max(s,swr),sor)-swr)/(sor-swr);
    fw = sm.^nw./(sm.^nw + (1-sm).^no);
    fp = sm.^nw./(sm.^nw + M*(1-sm).^no);
else
    [swr, sor, ads] = deal(0.25,0.75,0);
    % [nw,no,M] = deal(2,3,1);
    [nw,no,M] = deal(2,2.5,0.225);
    [swrs,sors,Msf] = deal(.05,.975,4);
    s  = linspace(0,1,501);
    sm = (min(max(s,swr),sor)-swr)/(sor-swr);
    fw = sm.^nw./(sm.^nw + M*(1-sm).^no);
    sm = (min(max(s,swrs),sors)-swrs)/(sors - swrs);
    fp = sm.^nw./(sm.^nw + Msf*M*1*(1-sm).^no);
    %load ff-sft2.mat;
    %[s, fw, fp] = deal(sW', fw', fs');
end
subplot(1,2,1),cla,hold all
plot(s,fw,'-','LineWidth',2);
plot(s,fp,'-','LineWidth',2);

%% Find envelope for polymer curve
i3 = find(fp(2:end)./(s(2:end)+ads)-diff(fp)./diff(s)>0,1);
plot([-ads s(i3:end)],fp([1 i3:end]),'--','LineWidth',2);
plot(s(i3),fp(i3),'o','LineWidth',2);

%% Find intersection of this envelope and the pure water curve
i2 = find(min(fp(i3)/(s(i3)+ads).*(s(2:end)+ads),1)<fw(2:end),1);
plot(s(i2),fw(i2),'o','LineWidth',2);

%% Find envelope for pure water curve
i0 = find(s<swr,1,'last')+1;
flag = fw(i0:i2)./(s(i0:i2)-swr)>diff(fw(i0:i2+1))./diff(s(i0:i2+1));
i1 = min(find(flag)+i0,i2);
if isempty(i1), i1=i2; end
plot(s(i1(1)),fw(i1(1)),'o','LineWidth',2)
plot(s([i0 i1:i2]),fw([i0 i1:i2]),'--','LineWidth',2)

%% Plot self-similar solution
subplot(1,2,2)
hold all
xi = fliplr(diff(fp(i3:end))./diff(s(i3:end)));
ss = fliplr(s(i3:end-1));
df = fliplr(diff(fw([i0 i1:i2]))./diff(s([i0 i1:i2])));
ds = fliplr(s(i1:i2));
xi = [xi([1:end end]), df];
ss = [ss, ds([1 1:end])];
xi = [xi([1:end end]), max(xi)+.5];
ss = [ss, s(i0), swr];
if polymer
    plot(xi(ss<=sor),ss(ss<=sor),'-','LineWidth',2);
else
    plot(xi,ss,'-','LineWidth',1);
end
axis tight, set(gca,'Ylim',[0 1]);

%% Plot solution without polymer
flag = fw(i0:end)./(s(i0:end)-swr)-eps> ...
    diff(fw([i0:end end]))./diff(s([i0:end end]));
iw  = min(find(flag,1,'first')+i0,numel(fw));
xiw = fliplr(diff(fw([i0 iw:end]))./diff(s([i0 iw:end])));
ssw = s(end:-1:iw);
xiw = [xiw xiw(end) max(xiw)+.5];
ssw = [ssw, swr swr];
plot(xiw(ssw<=sor),ssw(ssw<=sor),'--k','LineWidth',1);

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
