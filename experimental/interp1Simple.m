function varargout = interp1Simple(X, Y, xi)
%INTERP1 1-D interpolation (table lookup)
%
%   Some features of INTERP1 will be removed in a future release.
%   See the R2012a release notes for details.
%
%   Vq = INTERP1(X,V,Xq) interpolates to find Vq, the values of the
%   underlying function V=F(X) at the query points Xq. X must
%   be a vector of length N.
%   If V is a vector, then it must also have length N, and Vq is the
%   same size as Xq.  If V is an array of size [N,D1,D2,...,Dk], then
%   the interpolation is performed for each D1-by-D2-by-...-Dk value
%   in V(i,:,:,...,:).
%   If Xq is a vector of length M, then Vq has size [M,D1,D2,...,Dk].
%   If Xq is an array of size [M1,M2,...,Mj], then Vq is of size
%   [M1,M2,...,Mj,D1,D2,...,Dk].
%
%   Vq = INTERP1(V,Xq) assumes X = 1:N, where N is LENGTH(V)
%   for vector V or SIZE(V,1) for array V.
%
%   Interpolation is the same operation as "table lookup".  Described in
%   "table lookup" terms, the "table" is [X,V] and INTERP1 "looks-up"
%   the elements of Xq in X, and, based upon their location, returns
%   values Vq interpolated within the elements of V.
%
%   Vq = INTERP1(X,V,Xq,METHOD) specifies alternate methods.
%   The default is linear interpolation. Use an empty matrix [] to specify
%   the default. Available methods are:
%
%     'nearest'  - nearest neighbor interpolation
%     'linear'   - linear interpolation
%     'spline'   - piecewise cubic spline interpolation (SPLINE)
%     'pchip'    - shape-preserving piecewise cubic interpolation
%     'cubic'    - same as 'pchip'
%     'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate and uses 'spline' if X is not equally
%                  spaced.
%
%   Vq = INTERP1(X,V,Xq,METHOD,'extrap') uses the interpolation algorithm
%   specified by METHOD to perform extrapolation for elements of Xq outside
%   the interval spanned by X.
%
%   Vq = INTERP1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of the
%   interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
%   EXTRAPVAL.  The default extrapolation behavior with four input arguments
%   is 'extrap' for 'spline' and 'pchip' and EXTRAPVAL = NaN for the other
%   methods.
%
%   PP = INTERP1(X,V,METHOD,'pp') will use the interpolation algorithm specified
%   by METHOD to generate the ppform (piecewise polynomial form) of V. The
%   method may be any of the above METHOD except for 'v5cubic'. PP may then
%   be evaluated via PPVAL. PPVAL(PP,Xq) is the same as
%   INTERP1(X,V,Xq,METHOD,'extrap').
%
%   For example, generate a coarse sine curve and interpolate over a
%   finer abscissa:
%       X = 0:10; V = sin(X); Xq = 0:.25:10;
%       Vq = interp1(X,V,Xq); plot(X,V,'o',Xq,Vq)
%
%   For a multi-dimensional example, we construct a table of functional
%   values:
%       X = [1:10]'; V = [ X.^2, X.^3, X.^4 ];
%       Xq = [ 1.5, 1.75; 7.5, 7.75]; Vq = interp1(X,V,Xq);
%
%   creates 2-by-2 matrices of interpolated function values, one matrix for
%   each of the 3 functions. Vq will be of size 2-by-2-by-3.
%
%   Class support for inputs X, V, Xq, EXTRAPVAL:
%      float: double, single
%
%   See also INTERP1Q, INTERPFT, SPLINE, PCHIP, INTERP2, INTERP3, INTERPN, PPVAL.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 5.41.4.21 $  $Date: 2012/05/11 20:15:03 $
%
% Determine input arguments.
% Work backwards parsing from the end argument.

% Set up the defaults
narginchk(2,5);
method = [];
extrapolate = false;
extrapval = NaN;
extrapvalimposed=false;
ndataarg = nargin; % Number of X,V,Xq args. Init to nargin and reduce.

if ischar(varargin{end})
   ischarendm1 = ischar(varargin{end-1});
   if strcmp(varargin{end},'pp')
      if (nargin ~= 4)
         error(message('MATLAB:interp1:ppOutput'))
      end
      method = sanitycheckmethod(varargin{end-1});
      % X and V should be vectors of equal length
      X = varargin{1};
      V = varargin{2};
      if isvector(V)
         orig_size_v = size(V);
         V = V(:); % Reorient not considered a resize
      else
         orig_size_v = size(V);
         n = orig_size_v(1);
         ds = orig_size_v(2:end);
         prodDs = prod(ds);
         V = reshape(V,[n prodDs]);
      end
      sanitycheck(X,V);
      X = X(:);
      if isscalar(X)
         error(message('MATLAB:interp1:NotEnoughPts'))
      end
      if any(diff(X)<0)
         [X, idx] = sort(X);
         V = V(idx,:);
      end
      griddedInterpolant(X,V(:,1));  % Use this to sanity check the input.
      pp = ppinterp(X, V, orig_size_v, method);
      varargout{1} = pp;
      return
   elseif strcmp(varargin{end},'extrap')
      if (nargin ~= 4 && nargin ~= 5)
         error(message('MATLAB:interp1:nargin'));
      end
      if isempty(varargin{end-1})==false && ischarendm1==false
         error(message('MATLAB:interp1:ExtrapNoMethod'));
      end
      extrapolate = true;
      method = varargin{end-1};
      ndataarg = nargin-2;
   else
      if ischar(varargin{end-1})
         error(message('MATLAB:interp1:InvalidSpecPPExtrap'))
      end
      method = varargin{end};
      ndataarg = nargin-1;
   end
elseif isscalar(varargin{end}) && ischar(varargin{end-1})
   extrapval = varargin{end};
   extrapvalimposed=true;
   method = varargin{end-1};
   ndataarg = nargin-2;
elseif isscalar(varargin{end}) && isempty(varargin{end-1}) && (nargin == 4 || nargin == 5)
   % default method via []
   extrapval = varargin{end};
   extrapvalimposed = true;
   ndataarg = nargin-2;
elseif isempty(varargin{end})
   % This is potentially ambiguous, the assumed intent is case I
   % I)    X, V, []   Empty query
   % II)   V, [], [] Empty query and empty method,
   % III)  V, Xq, [] Empty method
   if nargin ~= 3
      ndataarg = nargin-1;
   end
end
if isempty(method)
   method = 'linear';
else
   if method(1) == '*'
      method(1) = [];
   end
   switch lower(method(1))
      case 'n'
         method = 'nearest';
      case 'l'
         method = 'linear';
      case 's'
         method = 'spline';
      case {'p', 'c'}
         method = 'pchip';
      case 'v'  % 'v5cubic'
         method = 'cubic';
         if(extrapolate)
            warning(message('MATLAB:interp1:NoExtrapForV5cubic'));
         end
         extrapolate = false;
      otherwise
         error(message('MATLAB:interp1:InvalidMethod'));
   end
end
% Set up X, V, and Xq and sanity check the data
% At this point we have two possible scenarios
% (X,V,Xq) or (V,Xq) and V may not be a vector
% if ndataarg ~= 2 or  ndataarg ~=3, error

if ndataarg == 2
   V = varargin{1};
   if isvector(V)
      orig_size_v = size(V);
      V = V(:); % Reorient not considered a resize
   else
      orig_size_v = size(V);
      n = orig_size_v(1);
      ds = orig_size_v(2:end);
      prodDs = prod(ds);
      V = reshape(V,[n prodDs]);
   end
   Xq = varargin{2};
   X =(1:size(V,1))';
elseif ndataarg == 3
   X = varargin{1};
   if ~isnumeric(X)
      error(message('MATLAB:interp1:Xnumeric'));
   end
   V = varargin{2};
   if isvector(V)
      orig_size_v = size(V);
      V = V(:); % Reorient not considered a resize
   else
      orig_size_v = size(V);
      n = orig_size_v(1);
      ds = orig_size_v(2:end);
      prodDs = prod(ds);
      V = reshape(V,[n prodDs]);
   end
   X = X(:);
   if any(diff(X)<0)
      [X, idx] = sort(X);
      V = V(idx,:);
   end
   Xq = varargin{3};
else
   error(message('MATLAB:interp1:nargin'));
end

if isscalar(X)
   if isempty(Xq)
      varargout{1} = zeros(size(Xq));
      return
   end
end

if (strcmpi(method,'pchip') || strcmp(method,'spline')) && any(find(isnan(V)))
   Vq = Interp1DStripNaN(X,V,Xq,method);
else
   if strcmpi(method,'linear') && extrapolate
      [X, V] = LinearShiftBounds(X, V, Xq);
   elseif strcmpi(method,'nearest') && extrapolate
      [X, V] = NearestShiftBounds(X, V, Xq);
   end
   Vq = Interp1D(X,V,Xq,method);
end


if extrapvalimposed
   % Impose the extrap val; this is independent of method
   extptids = Xq < X(1) | Xq > X(end);
   Vq(extptids,:) = extrapval;
end
if isvector(V)
   % V is a vector so size(Vq) == size(Xq)
   siz_vq = size(Xq);
else
   if isvector(Xq)
      % V is not a vector but Xq is. Batch evaluation.
      siz_vq = [length(Xq) orig_size_v(2:end)];
   else
      % Both V and Xq are non-vectors
      siz_vq = [size(Xq) orig_size_v(2:end)];
   end
end

% Reshape result, possibly to an ND array
varargout{1} = cast(reshape(Vq,siz_vq),superiorfloat(X,V,Xq));

end % INTERP1


function [X, V] = LinearShiftBounds(X, V, Xq)
% Run a sanity check before altering X
griddedInterpolant(X,V(:,1));
% Linearly shift lower bound
xqmin = min(Xq(:));
xqmax = max(Xq(:));
if xqmin < X(1)
   sloperat = (X(1) - xqmin)/(X(2) - X(1));
   V0 = V(1,:) - (V(2,:)-V(1,:)).*sloperat;
   X(1) = xqmin;
   V(1,:) = V0;
end
if xqmax > X(end)
   sloperat = (xqmax-X(end))/(X(end)-X(end-1));
   Vn1 = V(end,:) +  (V(end,:)-V(end-1,:)).*sloperat;
   X(end) = xqmax;
   V(end,:) = Vn1;
end
end


function [X, V] = NearestShiftBounds(X, V, Xq)
xqmin = min(Xq(:));
xqmax = max(Xq(:));
xtype = class(X);
vtype = class(V);
if xqmin < X(1) && xqmax > X(end)
   Xnew = zeros(numel(X)+2,1);
   Vnew = zeros(size(V,1)+2,size(V,2));
   Xnew = cast(Xnew,xtype);
   Vnew = cast(Vnew,vtype);
   Xnew(1,:) = xqmin;
   Vnew(1,:) = V(1,:);
   Xnew(2:(end-1),:) = X;
   Vnew(2:(end-1),:) = V;
   Xnew(end,:) = xqmax;
   Vnew(end,:) = V(end,:);
   X = Xnew;
   V = Vnew;
elseif xqmin < X(1)
   Xnew = zeros(numel(X)+1,1);
   Vnew = zeros(size(V,1)+1,size(V,2));
   Xnew = cast(Xnew,xtype);
   Vnew = cast(Vnew,vtype);
   Xnew(1,:) = xqmin;
   Vnew(1,:) = V(1,:);
   Xnew(2:end,:) = X;
   Vnew(2:end,:) = V;
   X = Xnew;
   V = Vnew;
elseif xqmax > X(end)
   Xnew = zeros(numel(X)+1,1);
   Vnew = zeros(size(V,1)+1,size(V,2));
   Xnew = cast(Xnew,xtype);
   Vnew = cast(Vnew,vtype);
   Xnew(1:(end-1),:) = X;
   Vnew(1:(end-1),:) = V;
   Xnew(end,:) = xqmax;
   Vnew(end,:) = V(end,:);
   X = Xnew;
   V = Vnew;
end
end


function Vq = Interp1D(X,V,Xq,method)

Xqcol = Xq(:);
num_vals = size(V,2);
if (num_vals>1)
   Xext = {cast(X,'double'),(1:num_vals)'};
   Xqext = {cast(Xqcol,class(Xext{1})),Xext{2:end}};
else
   Xext = {X};
   Xqext = {Xqcol};
end
if ~strcmpi(method,'pchip')
   F = griddedInterpolant(Xext,V,method);
end
if strcmpi(method,'pchip') || any(find(~isfinite(V)))
   F = griddedInterpolant(X,V(:,1),method);
   if strcmpi(method,'cubic') && strcmpi(F.Method,'spline') && any(find(isnan(V)))
      Vq = Interp1DStripNaN(X,V,Xq,'spline');
      return
   end
   Vq = zeros(numel(Xqcol),num_vals);
   Vq(:,1) = F(Xqcol);
   for iv = 2:num_vals
      F.Values = V(:,iv);
      Vq(:,iv) = F(Xqcol);
   end
   return;
end
Vq = F(Xqext);
end

function Vq = Interp1DStripNaN(X,V, Xq,method)

Xqcol = Xq(:);
num_value_sets = 1;
numXq = numel(Xqcol);
if ~isvector(V)
   num_value_sets = size(V,2);
end

% Allocate Vq
Vq = zeros(numXq,num_value_sets);
nans_stripped = false;
for i = 1:num_value_sets
   numvbefore = numel(V(:,i));
   [xi, vi] = stripnansforspline(X,V(:,i));
   numvafter = numel(vi);
   if numvbefore > numvafter
      nans_stripped = true;
   end
   F = griddedInterpolant(xi,vi,method);
   if isempty(Xq)
      Vq(:,i) = Xqcol;
   else
      Vq(:,i) = F(Xqcol);
   end
end
if nans_stripped
   warning(message('MATLAB:interp1:NaNstrip'));
end
end





%-------------------------------------------------------------------------%

function pp = ppinterp(X,V, orig_size_v)
%PPINTERP ppform interpretation.
n = size(V,1);
ds = 1;
prodDs = 1;
if ~isvector(V)
   ds = orig_size_v(2:end);
   prodDs = size(V,2);
end

breaks = X.';
page1 = (diff(V)./repmat(diff(X),[1, prodDs])).';
page2 = (reshape(V(1:end-1,:),[n-1, prodDs])).';
coefs = cat(3,page1,page2);
pp = mkpp(breaks,coefs,ds);

% Even if method is 'spline' or 'pchip', we still need to record that the
% input data V was oriented according to INTERP1's rules.
% Thus PPVAL will return Vq oriented according to INTERP1's rules and
% Vq = INTERP1(X,Y,Xq,METHOD) will be the same as
% Vq = PPVAL(INTERP1(X,Y,METHOD,'pp'),Xq)
pp.orient = 'first';

end % PPINTERP





