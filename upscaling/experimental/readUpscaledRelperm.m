function [T_w, T_o] = readUpscaledRelperm(do_plot, bc, grav)
%[T_w, T_o] = readUpscaledRelperm(do_plot, 'fixed', norm(gravity)~=0);
% show vl and cl
%

mydir = '/data/statoil_upscaling/newgeocell/ref_curves';
mydir = '/data/statoil_upscaling/newgeocell/ref_curves';
txt = 'geocell_4satnum_anisotropic-upscalerelperm';

%
%mydir = '/data/statoil_upscaling/newgeocell';
%txt = 'geocell_4satnum-upscalerelperm';

if nargin == 0
   do_plot = true;
   bc = 'periodic';
   grav = false;
end

if grav
   cl = 'gcl';
else
   cl = 'cl';

end
if strcmp(bc, 'fixed')
   bctxt = 'fbc';

   numcol = 5;
   marker = {'--*', '-*'};
else
   bctxt = 'pbc';

   numcol = 11;
   marker = {'--', '-'};
end

files = {[txt,'_o_', bctxt,'_', cl, '_imb.txt'], ...
   [txt,'_w_', bctxt,'_', cl, '_imb.txt']; ...
   [txt,'_o_', bctxt,'_vl_imb.txt'], ...
   [txt,'_w_', bctxt,'_vl_imb.txt']}; ...


%figure
%marker = {'--', '-'};
for i = 1:size(files, 1)

fid = fopen(fullfile(mydir, files{i, 1}));
t_o = readRelPermTable(fid, 1, numcol);
T_o{i} = t_o{1}(:, 2:5);
fclose(fid);

fid = fopen(fullfile(mydir, files{i, 2}));
t_w = readRelPermTable(fid, 1, numcol);
T_w{i} = t_w{1}(:, 2:5);
fclose(fid);

if do_plot
plot(t_o{1}(:, 2), t_o{1}(:, 3), ['k', marker{i}]); hold on;
plot(t_o{1}(:, 2), t_o{1}(:, 4), ['k', marker{i}]);
plot(t_o{1}(:, 2), t_o{1}(:, 5), ['k', marker{i}]);

legend('krx', 'kry', 'krz');

plot(t_w{1}(:, 2), t_w{1}(:, 3), ['k', marker{i}]);
plot(t_w{1}(:, 2), t_w{1}(:, 4), ['k', marker{i}]);
plot(t_w{1}(:, 2), t_w{1}(:, 5), ['k', marker{i}]);

xlabel('s_w');
ylabel('kr');
title('CL (--) and VL (-), periodic bc (*) ')
end
end


