function norne = NorneGeostat

% This function is used to specify all geostatistical parameters used to
% generate an initial ensemble. Mean values for permeability, porosity,
% net-to-gross, region multipliers, oil-water contacts, and fault
% multipliers, are loaded or copied from the Statoil benchmark case. 
%
% For more information we refer to the paper: 
%
% Lorentzen, R., Luo, X., Bhakta, T., Valestrand, R.: "History Matching 
% the Full Norne Field Model Using Seismic and Production Data", 
% SPE Journal, January 2019. DOI: https://doi.org/10.2118/194205-PA.
% 
% We also ask that the above paper is cited in publications aided by this
% code.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C): IRIS (International Research Institute of Stavanger), 2017. 
% Contact: Rolf.Lorentzen@iris.no


dim=[46,112,22]; % field dimension
ldim=dim(1)*dim(2); % layer dimension
norne.dim = dim;

% actnum
str = fileread('ACTNUM_0704.prop');
matchStr = regexp(str,'(?<=ACTNUM)(.*)(?=/)','match');
s = regexprep(matchStr{1},'\s+',' ');
act = str2num(s)'; %#ok<*ST2NM>
norne.actnum = act;

% porosity
str = fileread('PORO_0704.prop');
matchStr = regexp(str,'(?<=PORO)(.*)(?=/)','match');
s = regexprep(matchStr{1},'\s+',' ');
p = str2num(s)'; %#ok<*ST2NM>
p(act==0) = [];
for nr=1:dim(3)
    index=ldim*(nr-1)+1:1:ldim*nr;
    values=p(sum(act(1:index(1)-1))+1:sum(act(1:index(end))));   
    meanv(nr)=mean(values); %#ok<*AGROW>
    stdv(nr)=std(values);
end
norne.poroMean = p;
norne.poroLayerMean=meanv';
norne.poroLayerStd=stdv';
norne.poroStd=0.05;
norne.poroLB=0.1;
norne.poroUB=0.4;
norne.poroRange=26;

% permeability
str = fileread('PERM_0704.prop');
matchStr = regexp(str,'(?<=PERMX)(.*)(?=/)','match');
s = regexprep(matchStr{1},'\s+',' ');
k = str2num(s)';
k = log(k);
k(act==0) = [];
for nr=1:dim(3)
    index=ldim*(nr-1)+1:1:ldim*nr;
    values=k(sum(act(1:index(1)-1))+1:sum(act(1:index(end))));    
    meanv(nr)=mean(values);
    stdv(nr)=std(values);
end
norne.permxLogMean = k;
norne.permxLayerLnMean=meanv';
norne.permxLayerStd=stdv';
norne.permxStd=1;
norne.permxLB=0.1;
norne.permxUB=10;
norne.permxRange=26;

% correlation between layers
for nr=1:dim(3)-1
    index=ldim*(nr-1)+1:1:ldim*nr;
    index2=ldim*(nr)+1:1:ldim*(nr+1);
    actlayer1=act(index);    
    actlayer2=act(index2);
    active=actlayer1.*actlayer2;
    values1=[k(sum(act(1:index(1)-1))+1:sum(act(1:index(end)))) ;...
        p(sum(act(1:index(1)-1))+1:sum(act(1:index(end))))];
    values2=[k(sum(act(1:index2(1)-1))+1:sum(act(1:index2(end)))) ;...
        p(sum(act(1:index2(1)-1))+1:sum(act(1:index2(end))))];
    v1=[actlayer1;actlayer1];
    v1(v1==1) = values1;
    v2=[actlayer2;actlayer2];
    v2(v2==1) = values2;
    co=corrcoef(v1(active==1), v2(active==1));
    corrWithNextLayer(nr)=co(1,2);
end
norne.corrWithNextLayer=corrWithNextLayer';

% correlation between porosity and permeability
norne.poroPermxCorr=0.7; 

% net-to-gross
str = fileread('NTG_0704.prop');
matchStr = regexp(str,'(?<=NTG)(.*)(?=/)','match');
s = regexprep(matchStr{1},'\s+',' ');
ntg = str2num(s)';
ntg(act==0) = [];
for nr=1:dim(3)
    index=ldim*(nr-1)+1:1:ldim*nr;
    values=ntg(sum(act(1:index(1)-1))+1:sum(act(1:index(end))));    
    meanv(nr)=mean(values);
    stdv(nr)=std(values);
end
norne.ntgMean = ntg;
norne.ntgLayerMean=meanv';
norne.poroNtgCorr = 0.6;
norne.ntgStd=0.1;
norne.ntgLB=0.01;
norne.ntgUB=1;
norne.ntgRange=26;

% rel-perm end-point scaling
norne.krwMean=1.15;
norne.krwLB=0.8;
norne.krwUB=1.5;
norne.krgMean=0.9;
norne.krgLB=0.8;
norne.krgUB=1;

% oil-water contact
norne.owcMean=[2692.0; 2585.5; 2618.0; 2400.0; 2693.3];
norne.owcLB=norne.owcMean-10;
norne.owcUB=norne.owcMean+10;

% region multipliers
norne.multregtLogMean=log10([0.0008; 0.1;0.05]);
norne.multregtStd=0.5;
norne.multregtLB=-5;
norne.multregtUB=0;

% z-multipliers
norne.z1Mean=-2;
norne.z1Std=0.5;
norne.z8Mean=-1.3;
norne.z8Std=0.5;
norne.z11Mean=-2;
norne.z11Std=0.5;
norne.z12Mean=-2;
norne.z12Std=0.5;
norne.z15Mean=-2;
norne.z15Std=1;
norne.z18Mean=-2;
norne.z18Std=1;
norne.zLB=-4;
norne.zUB=0;
norne.multzRange = 26;

% fault multipliers
str = fileread('FAULTMULT_AUG-2006.INC');
matchStr = regexp(str,'''(.[^'']*)/','tokens');
v = zeros(53,1);
for I = 1:53
    v(I) = str2num(matchStr{I}{1});
end
norne.multfltLogMean=log10(v);
norne.multfltStd=0.5;
norne.multfltLB=-5;
norne.multfltUB=2;







