try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    % The ADFNE folder should be put within the shale module folder
    Globals;
    warning('on');
catch
    %% Print warning when ADFNE is not available
    warning('adfne:missing', 'ADFNE or Statistics & Machine Learning Toolbox is not available!');
    warning('off', 'adfne:missing');
end

rng(12345); % the random number ensures that the same realization of 
          % fracture network is generated every time the code is run
tol = 0.01;
physdim = [3300 1300 250]*ft;
celldim = [98 15 3];
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);

try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    
    set1 = Field(DFN('dim',3,'n',350,'dir',45,'ddir',-100,'minl',4,'mu',10, ...
         'maxl',30,'bbx',[tol,tol,tol,physdim(1)-tol, ...
         physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
         'shape','l','q',4),'Poly');
    set2 = Field(DFN('dim',3,'n',350,'dir',-45,'ddir',-100,'minl',4,'mu',10, ...
         'maxl',30,'bbx',[tol,tol,tol,physdim(1)-tol, ...
         physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
         'shape','l','q',4),'Poly');
catch
    %% Load generated fractures when ADFNE is not available
    load('StochasticFracs.mat');
end
fracSet = [set1;set2]; 
numFplanes = numel(fracSet); 
startIdx = 1;
for i=1:numFplanes
    idxGlobal = i;
    fracplanes(idxGlobal).points = fracSet{i}(1:end-1,:);
end
numHFplanes = numel(set1);
numNFplanes = numel(set2); 
wells=[];

plotfracSystem_NF(G,fracplanes,numHFplanes,wells,'label',false)
view(-5,60)

%% This commented code requires ADFNE and the Statistics & Machine Learning Toolbox

try      
    % Plot stereonet if ADFNE is available
    
    o = Orientation([set1;set2]);
    figure,
    Stereonet([o.Dip],[o.Dir],'density',true,'marker','*','ndip',6,...
        'ndir',24,'cmap',@jet,'color','y');
catch
    % Don’t plot if ADFNE is not available
end
