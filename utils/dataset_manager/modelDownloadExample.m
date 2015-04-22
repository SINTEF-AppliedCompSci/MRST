% Listing data directory
ddir = '/data/tmpdatasets';
!rm -rf /data/tmpdatasets
mkdir(ddir)
mrstDataDirectory(ddir)
% Setting data directory
% mrstDataDirectory('some/path')
addpath('gui')
addpath('utils')
addpath('datasets')
%% Get info struct for dataset
[info, present] = dataset_saigup();
info

% Download a dataset
mrstVerbose on
downloadDataset('saigup');

downloadDataset('norne');

%%
[info, present] = getAvailableDatasets();

%%
close all
mrstDatasetGUI()

%%
listDatasetExamples('spe1')

%%
getDatasetPath('Norne', 'download', false)

%%
downloadAllDatasets