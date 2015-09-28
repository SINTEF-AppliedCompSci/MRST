function [ plume ] = getLayer9CO2plumeOutlines( )

% the following .mat files which are loaded must be placed on the current
% working directory path. These files can be downloaded from 
% https://bitbucket.org/mrst/mrst-co2lab/downloads

% note that the years 2006a and 2006b correspond to the same plume data for
% the given year, however are separated into two different arrays in order
% to avoid the plotting of a connected line between separate polygons.

% years 1999, 2004, and 2008 also have more than one polygon to represent
% plume, thus two separate .mat files were created and are loaded here.


% 1999
load layer9_polygons_1999a.mat;
plume{1}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{1}.year    = 1999;

load layer9_polygons_1999b.mat;
plume{2}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{2}.year    = 1999;


% 2001
load layer9_polygons_2001.mat;
plume{3}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{3}.year    = 2001;


% 2002
load layer9_polygons_2002.mat;
plume{4}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{4}.year    = 2002;


% 2004
load layer9_polygons_2004a.mat;
plume{5}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{5}.year    = 2004;

load layer9_polygons_2004b.mat;
plume{6}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{6}.year    = 2004;


% 2006
load layer9_polygons_2006a.mat;
plume{7}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{7}.year    = 2006;

load layer9_polygons_2006b.mat;
plume{8}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{8}.year    = 2006;


% 2008
load layer9_polygons_2008a.mat;
plume{9}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{9}.year    = 2008;

load layer9_polygons_2008b.mat;
plume{10}.outline = CO2plumeOutline; clear CO2plumeOutline
plume{10}.year    = 2008;


end

