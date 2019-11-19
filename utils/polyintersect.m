function  varargout = polyintersect(varargin)
%POLYINTERSECT  Finds all intersections of 2 polygons.
%
%    [xr,yr] = polyintersect(x1,y1,x2,y2)
%
% returns all intersections between the line pieces of polygons
% (x1,y1) and (x2,y2).
%
% POLYINTERSECT calculates analytically any possible intersection
% if the line segments would have infinate length, and subsequently
% discards those intersections that are outside the
% support points (begin and end points) of each segment. This
% method is in principle fool-proof, (and also works for vertical line segments).
%
% Note that whether line 1 or line 2 is longest does
% not matter (nor for speed, not for solution).
% Note that line 1 and line2 can be NaN-separated polygons.
% Note that for two parallel lines no crossing is returned
% (even if they partly overlap).
%
%    [xr,yr] = polyintersect(x1,y1,x2,y2,<keyword,value>)
%
% implemented <keyword,value> pairs are:
%   * debug    :   0, 1 (default) for plot of intersects
%                  2 for plotting all segments of line 2 one by one.
%   * disp     :   0, 1 (default) for printing progress to command line
%
% See also: POLYJOIN, POLYSPLIT, (mapping toolbox)
%           POLY2STRUC,
%           FINDCROSSINGSOFPOLYGONANDPOLYGON, LANDBOUNDARY

%   --------------------------------------------------------------------
%   Copyright (C) 2008 Deltares
%       Gerben J. de Boer
%
%       gerben.deboer@deltares.nl
%
%       Deltares (former Delft Hydraulics)
%       Rotterdamseweg 185
%       2629 HD Delft
%       The Netherlands
%
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
%   USA
%   --------------------------------------------------------------------

%% For vertical lines where slope becomes 1/0=Inf,
%% ignire error messages, the Inf is
%% taken care off afterwards.
%% ----------------------

% tmp = warning('QUERY', 'MATLAB:divideByZero');
% if strcmp(tmp.state,'on')
%     warning off MATLAB:divideByZero
% end

%% Set defaults for keywords
%% ----------------------

OPT.debug = 0;
OPT.disp  = 0;

%% Return defaults
%% ----------------------

if nargin==0
    varargout = {OPT};
    return
else
    line(1).cor.x = varargin{1};
    line(1).cor.y = varargin{2};
    line(2).cor.x = varargin{3};
    line(2).cor.y = varargin{4};
end

%% Cycle keywords in input argument list to overwrite default values.
%% Align code lines as much as possible to allow for block editing in textpad.
%% Only start <keyword,value> pairs after the REQUIRED arguments.
%% ----------------------

iargin    = 5;
while iargin<=nargin,
    if      ischar(varargin{iargin}),
        switch lower(varargin{iargin})
            case 'debug'     ;iargin=iargin+1;OPT.debug    = varargin{iargin};
            case 'disp'      ;iargin=iargin+1;OPT.disp     = varargin{iargin};
            otherwise
                error(['Invalid string argument: ' varargin{iargin}]);
        end
    end;
    iargin=iargin+1;
end;

%% y = a*x + b
%% Get slope a and offset b
%% Note number of line pieces        (a,b) at line centers (# meshes  )
%% is one less then number of pionts (x,y) at line corners (# vertices)
%% ---------------------

for ii=1:length(line)
    line(ii).cen.a = diff(line(ii).cor.y)./....
        diff(line(ii).cor.x);        % is Inf for vertical lines
    line(ii).cen.b = line(ii).cor.y(1:end-1) - ...
        line(ii).cor.x(1:end-1).*...
        line(ii).cen.a;              % is NaN or -Inf for vertical lines
    line(ii).n_meshes = length(line(ii).cen.a);
end

%% Fill a matrix with analytical solution of all crossings
%% of (extended) line pieces (does not work if one line is vertical):
%%
%%           yc = a1*xc+b1
%%           yc = a2*xc+b2
%%
%%     a1*xc+b1 = a2*xc+b2
%%   (a1-a2)*xc =  b2-b1
%%           xc = (b2-b1)/(a1-a2)
%% ---------------------

%% preallocate with a certain arbitrary size
%% ---------------------

crosspoint.x = nan([max(line(1).n_meshes,line(2).n_meshes) 1]);
crosspoint.y = nan([max(line(1).n_meshes,line(2).n_meshes) 1]);
crosspoint.n = 0;

%-% NOTE: stroring antire matrix of (all pieces line1) x (all pieces line2) NOT POSSIBLE
%-% crosspoint.x = repmat(nan,[line(1).n_meshes line(2).n_meshes]); % this becomes too big with 2 long polygons
%-% crosspoint.y = repmat(nan,[line(1).n_meshes line(2).n_meshes]); % this becomes too big with 2 long polygons

for imesh1=1:line(1).n_meshes
    if OPT.disp
        disp(['Processing intersect line piece ',num2str(imesh1),' / ',num2str(line(1).n_meshes),' of line 1 with all pieces of line 2.'])
    end
    for imesh2=1:line(2).n_meshes
        
        local.x = (line(2).cen.b(imesh2) - line(1).cen.b(imesh1))./...
            (line(1).cen.a(imesh1) - line(2).cen.a(imesh2));
        local.y =  line(1).cen.a(imesh1).*local.x + ...
            line(1).cen.b(imesh1);
        
        %% Remove crossings outside line-piece:
        %%
        %%
        %%
        %%     @
        %%      \  crossing that comes from analytical solution, but that is outside end pionts of line 1, so discard to NaN!
        %%       \
        %%     ...#.........o------------------o line 1 (2 line pieces, 3 vertices)
        %%         \       /
        %%          \     /
        %%           \   /
        %%            \ /
        %%             @  crossing that comes OK from analytical solution
        %%            / \
        %%           /   \
        %%          /     \
        %%         /       \
        %%        /         \
        %%       /           \
        %%     o              \
        %%                     @ line 2 (1 line piece, 2 vertices)
        %%
        %% ---------------------
        
        %% Calculate crossings if one line segment is vertical
        %% ---------------------
        
        if isinf(local.x) || isnan(local.x)
            
            %% Deal with two vertical lines that do not yield a crossing (both they don't overlap and when they don't overlap)
            %% ---------------------
            
            if  isinf(line(1).cen.a(imesh1)) && ...
                    isinf(line(2).cen.a(imesh2))
                
                if OPT.disp==1
                    disp('2 parallel vertical lines encountered, any overlap dismissed.')
                end
                local.x = NaN;
                local.y = NaN;
                
                %% Deal with two horizontal lines that do not yield a crossing (both they don't overlap and when they don't overlap)
                %% ---------------------
            elseif  line(1).cen.a(imesh1)==line(2).cen.a(imesh2)
                if OPT.disp==1
                    disp('2 parallel lines encountered, any overlap dismissed.')
                end
                local.x = NaN;
                local.y = NaN;
                
                %% Deal with one vertical line
                %% ---------------------
            elseif  isinf(line(1).cen.a(imesh1))
                if OPT.disp==1
                    disp('vertical lines encountered polygon 1.')
                end
                local.x = line(1).cor.x(imesh1);              % get x from vertical line
                local.y = line(2).cen.a(imesh2)*local.x + ... % get y from slope other line
                    line(2).cen.b(imesh2);
                
                %% Deal with other vertical line
                %% ---------------------
            elseif isinf(line(2).cen.a(imesh2))
                if OPT.disp==1
                    disp('vertical lines encountered polygon 2.')
                end
                local.x = line(2).cor.x(imesh2);              % get x from vertical line
                local.y = line(1).cen.a(imesh1)*local.x + ... % get y from slope other line
                    line(1).cen.b(imesh1);
            end
            
        end % if isinf(local.x) | isnan(local.x)
        
        %% Check whether crossing is on line segments
        %% Deal with with all other lines (incl. two parallel ones)
        %% ---------------------
        
        %NEWCODE start (probably faster)
        eps=1e-4;
        if        local.x < min(line(1).cor.x(imesh1 + [0 1]))-eps; %   *   + ----+
            local.x = NaN;
            local.y = NaN;
        else %1
            if       local.x > max(line(1).cor.x(imesh1 + [0 1]))+eps; %       + ----+   *
                local.x = NaN;
                local.y = NaN;
            else %2
                if      local.y < min(line(1).cor.y(imesh1 + [0 1]))-eps; %   *   + ----+      @ 90 DEG
                    local.x = NaN;
                    local.y = NaN;
                else %3
                    if     local.y > max(line(1).cor.y(imesh1 + [0 1]))+eps; %       + ----+   *  @ 90 DEG
                        local.x = NaN;
                        local.y = NaN;
                    else %4
                        if    local.x < min(line(2).cor.x(imesh2 + [0 1]))-eps; %   *   + ----+
                            local.x = NaN;
                            local.y = NaN;
                        else %5
                            if   local.x > max(line(2).cor.x(imesh2 + [0 1]))+eps; %       + ----+   *
                                local.x = NaN;
                                local.y = NaN;
                            else %6
                                if  local.y < min(line(2).cor.y(imesh2 + [0 1]))-eps; %   *   + ----+      @ 90 DEG
                                    local.x = NaN;
                                    local.y = NaN;
                                else %7
                                    if local.y > max(line(2).cor.y(imesh2 + [0 1]))+eps; %       + ----+   *  @ 90 DEG
                                        local.x = NaN;
                                        local.y = NaN;
                                    else %8
                                        
                                        if ~isnan(local.x) && OPT.debug>1
                                            disp(['found crosssection in polygon 1 section ',num2str(imesh1),' polygon 2 section ',num2str(imesh2)])
                                        end
                                        
                                    end %8
                                end %7
                            end %6
                        end %5
                    end %4
                end %3
            end %2
        end %1
        
        %NEWCODE end (probably faster)
        
        %OLDCODE start (probably slower)
        %            if  local.x < min(line(1).cor.x(imesh1 + [0 1])) | ... %   *   + ----+
        %                local.x > max(line(1).cor.x(imesh1 + [0 1])) | ... %       + ----+   *
        %                local.y < min(line(1).cor.y(imesh1 + [0 1])) | ... %   *   + ----+      @ 90 DEG
        %                local.y > max(line(1).cor.y(imesh1 + [0 1])) | ... %       + ----+   *  @ 90 DEG
        %                local.x < min(line(2).cor.x(imesh2 + [0 1])) | ... %   *   + ----+
        %                local.x > max(line(2).cor.x(imesh2 + [0 1])) | ... %       + ----+   *
        %                local.y < min(line(2).cor.y(imesh2 + [0 1])) | ... %   *   + ----+      @ 90 DEG
        %                local.y > max(line(2).cor.y(imesh2 + [0 1]));      %       + ----+   *  @ 90 DEG
        %
        %                local.x = NaN;
        %                local.y = NaN;
        %
        %              %else
        %              %   if ~isnan(local.x)
        %              %      disp(['found ',num2str(imesh2)])
        %              %   end
        %              %end
        %
        %            end
        %OLDCODE end (probably slower)
        
        %% Keep only real crossings
        %% ---------------------
        mask         = ~isnan(local.x) & ...
            ~isnan(local.y);
        local.x = local.x(mask);
        local.y = local.y(mask);
        local.n = length(local.x);
        
        %% Store local real crossings into big vector
        %% ---------------------
        
        index1                      = crosspoint.n + 1;
        index2                      = crosspoint.n + local.n;
        crosspoint.x(index1:index2) = local.x;
        crosspoint.y(index1:index2) = local.y;
        crosspoint.n                = crosspoint.n + local.n;
        
        %% Plot local real crossings for debug
        %% ---------------------
        
        if OPT.debug > 1
            TMP = figure;
            plot(line(2).cor.x,...
                line(2).cor.y,'b.-')
            axis equal
            %axis(axis)
            hold on
            grid on
            plot(line(1).cor.x,...
                line(1).cor.y,'k-')
            
            plot(line(2).cor.x(imesh2 + [0 1]),...
                line(2).cor.y(imesh2 + [0 1]),'b.-','linewidth',5)
            plot(line(1).cor.x(imesh1 + [0 1]),...
                line(1).cor.y(imesh1 + [0 1]),'k.-','linewidth',5)
            
            plot(local.x,...
                local.y,'ro')
            plot(local.x,...
                local.y,'r.')
            title(['POLYINTERSECT: All crosssings of 1st polygon ',num2str(imesh1),', with 2nd polygon line piece: ',num2str(imesh2)])
            disp('Press key to continue ...')
            pause
            try
                close(TMP)
            catch
                disp('');
            end
        end % debug
        
    end % imesh2
    
end % imesh1

%% Keep only real crossings (remove ones that were allocated to much)
%% ---------------------

mask         = ~isnan(crosspoint.x) & ...
    ~isnan(crosspoint.x);
crosspoint.x = crosspoint.x(mask);
crosspoint.y = crosspoint.y(mask);

%% Debug (optional)
%% ---------------------

if OPT.debug > 0
    TMP = figure;
    plot(line(2).cor.x,...
        line(2).cor.y,'b.-')
    
    axis equal
    axis(axis)
    hold on
    grid on
    plot(line(1).cor.x,...
        line(1).cor.y,'k.-')
    plot(crosspoint.x,...
        crosspoint.y,'ro')
    plot(crosspoint.x,...
        crosspoint.y,'r.')
    title('POLYINTERSECT: All crosssings of entire 1st polygon with entire 2nd polygon.')
    disp('Press key to continue ...')
    pause
    try
        close(TMP)
    catch
        disp('');
    end
end % debug

%% Output
%% ---------------------

if nargout==2
    varargout = {crosspoint.x,crosspoint.y};
elseif nargout==3
    varargout = {crosspoint.x,crosspoint.y,OPT};
end

%% Restore previous state
%% ----------------------

% if strcmp(tmp.state,'on')
%     warning on MATLAB:divideByZero % put back on (it if were), after switching it off (if it weren't)
% end

%% EOF