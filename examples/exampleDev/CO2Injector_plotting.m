close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui
mrstVerbose on

open('CO2Injector_run.mat')

nx = 100;
nz = 20;

x = 1:nx;
z = 1:nz;





v = VideoWriter('pHcontour.avi');
open(v);

figure(1); box on;
xlabel('x')
ylabel('y');
set(gca,'nextplot','replacechildren'); 
drawnow;

for i = 1 :5: numel(states)
phmat = -log10(reshape(states{i}.components(:,1), nx,nz));

contourf(x,z, phmat');  
   colorbar;
   set(gca, 'ydir', 'reverse')
frame = getframe(gcf);
   writeVideo(v,frame);

end

close(v);
