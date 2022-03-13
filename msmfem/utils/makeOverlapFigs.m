%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
require coarsegrid 
printpath = fullfile(fileparts(mfilename('fullpath')), 'html');

f = figure;
set(gcf,'PaperPositionMode','auto');
pos = get(gcf,'position'); pos(3:4)=0.3*pos(3:4); set(gcf,'position',pos);
%%
%
%                  . . . . . . . . . .
%                . o-----------------o .
%                . |        |        | .
%                . |   i    |   j    | .
%                . |        |        | .
%                . o-----------------o .
%                  . . . . . . . . . .
[x,y]=meshgrid(-10:10,-9:9);
[Xc,Yc]=meshgrid(-10:5:10,-9:6:9);
newplot, hold on
   patch([-6 6 6 -6],[-3 -3 3 3],[0.85 0.85 0.85]);
   patch([-5 5 5 -5],[-4 -4 4 4],[0.85 0.85 0.85]);
   plot(x,y,'k',x',y','k');
   plot(Xc,Yc,'k',Xc',Yc','k','LineWidth',2);
   text(-2.5, 0, 'i','FontSize',32);
   text( 2.5, 0, 'j','FontSize',32);
hold off
axis([-8 8 -6 6]), axis off
if isdir(printpath),
   print('-dpng', fullfile(printpath, 'basisIJ'));
end
%%
%                           . . . . .
%                           .       .
%                      o--------o   .
%                      |       *|   .
%                      |        | . .
%                      |        |
%                      o--------o
%
newplot, hold on
   patch([-5 0 0 -5], [-3 -3 3 3],[0.85 0.85 0.85]);
   %  level 1
   patch([-1 0 0 -1], [3 3 4 4],[0.75 0.875 0.75]); text(-0.75,3.5,'1','FontSize',6);
   patch([ 0 1 1  0], [2 2 3 3],[0.75 0.875 0.75]); text( 0.25,2.5,'1','FontSize',6);
   %  level 2
   patch([-2 -1 -1 -2], [3 3 4 4],[0.875 0.75 0.75]); text(-1.75,3.5,'2','FontSize',6);
   patch([ 0  1  1  0], [1 1 2 2],[0.875 0.75 0.75]); text( 0.25,1.5,'2','FontSize',6);
   patch([-1  0  0 -1], [4 4 5 5],[0.875 0.75 0.75]); text(-0.75,4.5,'2','FontSize',6);
   patch([ 1  2  2  1], [2 2 3 3],[0.875 0.75 0.75]); text( 1.25,2.5,'2','FontSize',6);
   patch([ 0  1  1  0], [3 3 4 4],[0.875 0.75 0.75]); text( 0.25,3.5,'2','FontSize',6);
   %  level 3
   patch([-3 -2 -2 -3], [3 3 4 4],[0.75 0.75 0.875]); text(-2.75,3.5,'3','FontSize',6);
   patch([ 0  1  1  0], [0 0 1 1],[0.75 0.75 0.875]); text( 0.25,0.5,'3','FontSize',6);
   patch([-2 -1 -1 -2], [4 4 5 5],[0.75 0.75 0.875]); text(-1.75,4.5,'3','FontSize',6);
   patch([ 1  2  2  1], [1 1 2 2],[0.75 0.75 0.875]); text( 1.25,1.5,'3','FontSize',6);
   patch([ 1  2  2  1], [3 3 4 4],[0.75 0.75 0.875]); text( 1.25,3.5,'3','FontSize',6);
   patch([ 0  1  1  0], [4 4 5 5],[0.75 0.75 0.875]); text( 0.25,4.5,'3','FontSize',6);
   patch([ 2  3  3  2], [2 2 3 3],[0.75 0.75 0.875]); text( 2.25,2.5,'3','FontSize',6);
   patch([-1  0  0 -1], [5 5 6 6],[0.75 0.75 0.875]); text(-0.75,5.5,'3','FontSize',6);
  %
   plot(x,y,'k',x',y','k');
   plot(Xc,Yc,'k',Xc',Yc','k','LineWidth',2);
   plot(-0.5, 2.5, 'r.', 'MarkerSize',12);
hold off
axis([-7 4 -5 7]), axis off
if isdir(printpath),
   print('-dpng', fullfile(printpath, 'basisOverlapWell'));
end
%%                     . . . . . .
%                    . o---------o .
%                    . |         | .
%                    . |    *    | .
%                    . |         | .
%                    . o---------o .
%                      . . . . . .
newplot, hold on
   patch([-5 0 0 -5], [-4 -4 4 4],[0.85 0.85 0.85]);
   patch([-6 1 1 -6], [-3 -3 3 3],[0.85 0.85 0.85]);
   plot(x,y,'k',x',y','k');
   plot(Xc,Yc,'k',Xc',Yc','k','LineWidth',2);
   plot(-0.5, 2.5, 'r.', 'MarkerSize',12);
hold off
axis([-7 4 -5 7]), axis off
if isdir(printpath),
   print('-dpng', fullfile(printpath, 'basisOverlapBlock'));
end
%%
close(f);
