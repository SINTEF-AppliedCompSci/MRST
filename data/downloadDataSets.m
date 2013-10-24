function downloadDataSets(varargin)
%% Download data sets for CO2 laboratory
% This script will download and install all data sets used by the CO2
% laboratory that are publicly available and do not require registration by
% the user. The following data sets will be downloaded:
% 
% * The Johansen formation is a candidate site for large-scale CO2 storage
%   offshore the south-west coast of Norway. The Johansen data set,
%   consists of four different sector models and a full-field model, see
%   the <http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/
%   MatMoRA webpage> for more details. Here, we only download the NPD5
%   sector model.
%
% * The IGEMS conceptual model, used in one of the VE tutorials, is a
%   syntetic model set that was created early in the
%   <http://www.nr.no/igems IGEMS project> as an example of a large 3D grid
%   model with variations in the topsurface topography.
%
% In addition, a few examples use the Sleipner data set from ieaghg.org,
% which requires user registration and hence must be downloaded manually;
% see "sleipner/README" for more details.

if nargin > 0
    ask = false;
    given = @(x) any(strcmpi(varargin{1}, x));
    if given('all')
        getJohansen = true;
        getSloping = true;
        getAtlas = true;
        getIGEMS = true;
    else
        getJohansen = given('johansen');
        getSloping = given('slopingaquifer');
        getAtlas = given('atlas');
        getIGEMS = given('igems');
    end
else
    ask = true;
    getJohansen = false;
    getSloping = false;
    getAtlas = false;
    getIGEMS = false;
end

% Johansen model
if getJohansen || (ask && userConsent('Download the Johansen data set (~13 MB)'))
   fprintf(1, '  Downloading Johansen sector model (NPD5) ..');
   unzip('http://www.sintef.no/project/MatMoRA/Johansen/NPD5.zip', ...
      fullfile(VEROOTDIR, 'data', 'johansen'));
   fprintf(1,'done\n');

   fprintf(1,'  Downloading Johansen full-field model..');
   unzip('http://www.sintef.no/project/MatMoRA/Johansen/FULLFIELD_Eclipse.zip', ...
      fullfile(VEROOTDIR, 'data', 'johansen'));
    fprintf(1,'done\n');  
end

% IGEMS conceptual model
if getSloping || (ask && userConsent('Download synthetic sloping aquifer model (~3 MB)'))
   fprintf(1,'  Downloading synthetic sloping aquifer model ..');
   untar('http://www.sintef.no/project/MRST/IGEMS.tar.gz',...
      fullfile(VEROOTDIR, 'data'));
   fprintf(1,'done\n');
end

% CO2 Atlas
if getAtlas || (ask && userConsent('Download CO2 Atlas (~14 MB)'))
   fprintf(1,'  Downloading CO2 Atlas (this may take a few minutes) ..');
   untar('http://www.sintef.no/project/MRST/CO2_atlas.tar.gz',...
      fullfile(VEROOTDIR, 'data'));
   untar('http://www.sintef.no/project/MRST/mapAndWells.tar.gz', ...
      fullfile(VEROOTDIR, 'data', 'atlas'));
   fprintf(1,'done\n');
end

% IGEMS data set
if getIGEMS || (ask && userConsent('Download IGEMS data set (~1.1 GB)'))
   % Download and check surfaces
   fprintf(1,'  Downloading IGEMS surfaces (this may take several minutes) ..');
   datadir = fullfile(VEROOTDIR, 'data', 'igems');
   unzip('http://www.nr.no/files/sand/Igems/surfaces.zip', datadir);
   dirs = dir(fullfile(datadir,'surfaces'));
   for i=1:numel(dirs)
      if sum(isletter(dirs(i).name))==0, continue, end
      if numel(dir(fullfile(datadir,'surfaces',dirs(i).name,'*.irap')))<100
         disp(['Installation of surfaces/' dirs(i).name ' seems wrong']);
      end
   end
   fprintf(1,'done\n');
   
   % Download and check eclipse data sets
   fprintf(1, '  Downloading IGEMS grids (this may take several minutes) ..');
   unzip('http://www.nr.no/files/sand/Igems/eclipsegrids.zip', datadir);
   if numel(dir(fullfile(datadir,'eclipsegrids','*.GRDECL')))~=15
      disp('Installation of eclipsegrids/ seems wrong');
   end
   fprintf(1,'done\n');
end
