function appDir = AddMRST(mrstVersion)  
% <keywords>
%
% Purpose : Start the MRST package and add the required modules
%
% Syntax : 
%   appDir = AddMRST(mrstVersion);
%
% Input Parameters :
%   mrstVersion: string version of the MRST package to use e.g. 2020a
%
% Return Parameters :
%   The directory at which this functions exists to be used in the
%   configure functions
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
    % Change the current folder to the folder of this m-file.
    appDir = fileparts(which(mfilename));
    if(~isdeployed)
        cd(appDir);
        % SCAL software directory
        addpath(appDir);
        addpath(genpath(appDir + "\Source"));
        % MRST software directory
        mrstFolder = strcat("mrst-",mrstVersion);
        mrstDir = "../" + mrstFolder;
        % addpath(genpath(char(mrstDir)),"-frozen");
        addpath(mrstDir)
        % MRST startup routine
        mrststartup;
        mrstModule add incomp ad-core ad-blackoil ad-props mrst-gui
%         disp('----------------------------------------------')
%         disp('  MRST modules have been loaded successfully  ')
%         disp('----------------------------------------------')
    end
end