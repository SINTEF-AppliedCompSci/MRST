function generateNorneEnsemble(ensembleSize, randomNumber)

% This function generate an initial ensemble for the full Norne field.
% The ensemble is generated for porosity, log-permeability, net-to-gross,
% z-multipliers for base of layers [1,8,11,12,15,18], fault multipliers,
% oil-water contacts, multipliers for scaling of rel-perm curves, and
% region multipliers. Total number of generated parameters is 148159.
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
% Input:
% - ensembleSize : Number of generated ensemble members (defalult 100)
% - randomNumber : Number used to initialize random generator 
%                  (default 1.2345e5)
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


% set random generator
if nargin < 2
    rng(1.2345e5);
else
    if isscalar(randomNumber) && isreal(randomNumber) && ...
       rem(randomNumber,1) == 0 && randomNumber >= 0
        rng(randomNumber);
    else
        error('Random seed must be a nonnegative integer')
    end
end

% number of ensemble members
if nargin < 1
    N = 100;
else
    if isscalar(ensembleSize) && isreal(ensembleSize) && ...
       rem(ensembleSize,1) == 0 && ensembleSize > 0
        N = ensembleSize;
    else
        error('Ensemble size must be a positive integer')
    end
end
    
% set geostatistical parameters and initialize
norne = NorneGeostat();
A = norne.actnum; % mask for active gridcells
D = norne.dim; % reservoir dimention
N_L = D(1)*D(2); % number of gridcells in a layer
N_F = prod(D); % total number of gridcells

% mean values for poro, log-permx and ntg
M(:,1) = norne.poroMean; 
M(:,2) = norne.permxLogMean;
M(:,3) = norne.ntgMean;

% std for poro, permx and ntg
S(:,1) = norne.poroStd;
S(:,2) = norne.permxStd;
S(:,3) = norne.ntgStd;

% mean and std for multz
A_L = reshape(A,prod(D(1:2)),D(3)); % help variable (layerwise actnum)
M_M = [norne.z1Mean,norne.z8Mean,norne.z11Mean,norne.z12Mean,norne.z15Mean,norne.z18Mean];
S_M = [norne.z1Std,norne.z8Std,norne.z11Std,norne.z12Std,norne.z15Std,norne.z18Std];

% mean and std for multflt and multreg
M_MF = norne.multfltLogMean;
S_MF = norne.multfltStd;
M_MR = norne.multregtLogMean;
S_MR = norne.multregtStd;

% bounds for krw+krg
KR_LU = [norne.krwLB, norne.krwUB; norne.krgLB, norne.krgUB];

% bounds for owc
OWC_LU = [norne.owcLB,norne.owcUB];

% mean correlation lengths (ranges)
C(:,1) = norne.poroRange;
C(:,2) = norne.permxRange;
C(:,3) = norne.ntgRange;
if isfield('norne','multzRange')
    C_MZ = norne.multzRange;
else
    C_MZ = 26; 
end

% std for correlation lengths
if isfield('norne','C_S')
    C_S = norne.C_S;
else
    C_S = 2; 
end

% Poro / Permx correlation
R1 = norne.poroPermxCorr;

% Poro / NTG correlation
R2 = norne.poroNtgCorr;

% Lower/upper bound on relperm
norne.krwLB = [KR_LU(1,1)*ones(4,1); KR_LU(2,1)*ones(4,1)];
norne.krwUB = [KR_LU(1,2)*ones(4,1); KR_LU(2,2)*ones(4,1)];

% initial ensemble for Poro, Permx, NTG and Free parameters
ensemble(N) = struct('multz',[],'multflt',[],'multreg',[], ...
    'krw',[],'owc',[],'PORO',[],'PERMX',[],'NTG',[]);
for I = 1:N
    
    % multz
    A_MZ = A_L(:,[1 8 11 12 15 18]);
    A_MZ = A_MZ(:);
    V = ones(N_L,1);
    M_MZ = kron(M_M',V);
    S_MZ = kron(S_M',V);
    X = GaussianWithVariableParameters([D(1),D(2),6],M_MZ,S_MZ,C_MZ,C_S);
    ensemble(I).multz = min(max(X(A_MZ==1),norne.zLB),norne.zUB);
   
    % multflt
    N_MF = length(M_MF);
    X = M_MF + S_MF*randn(N_MF,1);
    ensemble(I).multflt = min(max(X,norne.multfltLB),norne.multfltUB);
    
    % multreg
    N_MR = length(M_MR);
    X = M_MR + S_MR*randn(N_MR,1);
    ensemble(I).multreg = min(max(X,norne.multregtLB),norne.multregtUB);
    
    % krw+krw
    X = [KR_LU(1,1) + (KR_LU(1,2)-KR_LU(1,1))*rand(4,1);...
         KR_LU(2,1) + (KR_LU(2,2)-KR_LU(2,1))*rand(4,1)];
    ensemble(I).krw = min(max(X,norne.krwLB),norne.krwUB);
    
    % owc
    N_OWC = 5;
    X = OWC_LU(:,1) + (OWC_LU(:,2)-OWC_LU(:,1)).*rand(N_OWC,1);
    ensemble(I).owc = min(max(X,norne.owcLB),norne.owcUB);
    
    % poro
    X1 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,1),C_S);
    ensemble(I).PORO = min(max(M(:,1) + S(:,1)*X1(A==1),norne.poroLB),norne.poroUB);
    
    % permx
    X2 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,2),C_S);
    X = R1*X1 + sqrt(1-R1^2)*X2;
    ensemble(I).PERMX = min(max(M(:,2) + S(:,2)*X(A==1),norne.permxLB),norne.permxUB);
    
    % ntg 
    X2 = GaussianWithVariableParameters(D,zeros(N_F,1),1,C(:,3),C_S);
    X = R2*X1 + sqrt(1-R2^2)*X2;
    ensemble(I).NTG = min(max(M(:,3) + S(:,3)*X(A==1),norne.ntgLB),norne.ntgUB);
   
end

% Save the ensemble

save(fullfile(getDatasetPath('norne_ensemble'),'data','NorneInitEns.mat'), 'ensemble');


%-----------------------------------------------
% Subfunction 
%----------------------------------------------- 
function [x,corrLength] = GaussianWithVariableParameters(fieldDim,meanValue, ...
			  Sdev,meanCorrLength,stdCorrLength)

% Setup a gaussian field with correlation length drawn from a
% normal distribution. The horizontal layers are genereated
% independently
%
% Input:
% - fieldDim       : dimension of the field.
% - meanValue      : the mean value of the field (vector of the size
%                    of the field).
% - Sdev           : standard deviation of the field.
% - meanCorrLength : mean correlation length.
% - stdCorrLength  : standard deviation of the correlation length.
%
% Output:
% - x              : the generated field.
% - corrLenght     : the drawn correlation length.

corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
if length(fieldDim)<3
    
    x=meanValue+FastGaussian(fieldDim,Sdev,corrLength);
    
else
    
    layerDim=prod(fieldDim(1:2));
    
    % Initialization
    x=meanValue;
    if length(Sdev)==1
        for I=1:fieldDim(3)
            x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
                I*layerDim)+FastGaussian(fieldDim(1:2),Sdev,corrLength);
            
            % Generate new correlation length for the next layer
            corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
        end
    else
        for I=1:fieldDim(3)
            x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
                I*layerDim)+ ...
                FastGaussian(fieldDim(1:2),Sdev(1+(I-1)*layerDim:I* ...
                layerDim),corrLength);
            
            % Generate new correlation length for the next layer
            corrLength=AddGnoise(meanCorrLength,stdCorrLength,1);
        end
    end
    
end



%-----------------------------------------------
% Subfunction 
%----------------------------------------------- 
function [Y,RTSIGMA] = AddGnoise( Ytrue,SIGMA,SQ )

% Add noise, normally distributed, with covariance given by SIGMA.
%
% Input:
% - Ytrue    Original signal.
% - SIGMA    Specification of covariance matrix of noise.  May be
%            entered as scalar, vector or full matrix. If it SIGMA
%            is a vector, then it is interpreted as the covariance
%            matrix is diag(SIGMA).
% - SQ       If present, determine whether SIGMA or SIGMA*SIGMA' is used
%            as the covariance matrix.  Thus, if the square root of
%            the covariance matrix has allready been calculated
%            previously, work may be saved by setting SQ to 1.
%
% Output:
% - Y        Signal with noise added.
% - RTSIGMA  The square root of SIGMA; RTSIGMA*RTSIGMA' = SIGMA.
%            (Helpful if it is cumbersome to compute).

% Compute the normally distributed noise, with covariance matrix
% SIGMA or SIGMA*SIGMA'.
try
    
    if nargin > 2 && SQ == 1
        % Use SIGMA*SIGMA' as covariance matrix
        RTSIGMA = SIGMA ;
        if min(size(SIGMA)) == 1
            % SIGMA is a scalar or vector
            error = RTSIGMA.*randn(size(Ytrue)) ;
        else
            error = RTSIGMA*randn(size(RTSIGMA,2),1) ;
        end
    else
        % Use SIGMA as covariance matrix
        if min(size(SIGMA))==1
            % SIGMA is entered as a scalar or a vector
            RTSIGMA = realsqrt(SIGMA);
            error = RTSIGMA.*randn(size(Ytrue));
        else
            [RTSIGMA,p] = chol(SIGMA); % The matrix must be transposed.
            if p>0
                disp('Problem with Cholesky factorization')
                disp(['p = ',num2str(p)]);
                RTSIGMA = real(sqrtm(SIGMA));
                disp('Finnaly - we got a square root!')
            end
            RTSIGMA = RTSIGMA';
            error = RTSIGMA*randn(size(Ytrue));
        end %if
    end
    % Add the noise:
    Y = Ytrue+error;
  
catch err
    
  disp('Error in AddGnoise');
  disp('Size Ytrue:')
  size(Ytrue)
  disp('Size SIGMA:')
  size(SIGMA)
  rethrow(err);
  
end

