function SaveVideo(frames,fileName)
%
% DESCRIPTION: saves a video from the simulation
%
% SYNOPSIS:
%   SaveResults(model)
%
% PARAMETERS:
%   - model - struct after the completion of the simulation
%   - fileName - name or path of the file to save the video
%
% RETURNS:
%   outputs a video from the simulation
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
axis tight manual 
set(gca,'nextplot','replacechildren');
format = 'Uncompressed AVI';
v = VideoWriter(fileName,format);
compressionRatio = 5; 
frameRate = 5;
quality = 100;
v.FrameRate = frameRate;
% v.CompressionRatio = compressionRatio;
% v.Quality = quality;
open(v)
writeVideo(v,frames);
close(v)