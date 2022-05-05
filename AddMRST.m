function AddMRST(mrstVersion)  
%
% DESCRIPTION: Add the required modules and starts the main MRST package
%
% SYNOPSIS:
%   AddMRST(mrstVersion);
%
% PARAMETERS:
%   mrstVersion - string containing MRST version
%
% EXAMPLE:
%   AddMRST("2020a") 
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
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