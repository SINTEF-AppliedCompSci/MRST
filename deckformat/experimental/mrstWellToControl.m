function control=mrstWellToControl(W,g,control,RUNSPEC,varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt=struct('add_wellindex', false);
opt = merge_options(opt,varargin{:});
% prototypes used
WELLSPECS_org={'INJE','I',[1],[1],[NaN],'WATER',[0],'STD','SHUT','NO',...
    [0],'SEG',[0];
    'PROD','P',[1],[1],[NaN],'OIL',[0],'STD','SHUT','NO',...
    [0],'SEG',[0]};

COMPDAT_org={'INJE01'    [1]    [1]    [1]    [1]    'OPEN'    [-1]...
    [-1]    [0.1000]    [-1]    [0]    'Default'    'Z'    [-1];...
    'PROD01'    [3]    [1]    [1]    [1]    'OPEN'    [-1]    [-1]... 
    [0.1000]    [-1]    [0]    'Default'    'Z'    [-1]
    };

v=Inf
WCONINJE_org={'Injector2'    'WATER'    'OPEN'    'rate'    [v] ...
    [v]    [v]    [v]    [0]    [v]    [v]    [v]    [v]    [v]};

WCONPROD_org={'Producer1'    'OPEN'    'bhp'    [v]    [v]    [v]    [v]...
    [v]    [v]    [v]    [0]    [v]};


cartDims=g.cartDims;
indexMap=g.cells.indexMap;
clear g

WELLSPECS=cell(numel(W),13);
for i=1:numel(W)
    [xi,xj,xk]=ind2sub(cartDims, indexMap(W(i).cells));
    for k=1:13
        if(W(i).sign>0)
                WELLSPECS{i,k} = WELLSPECS_org{1,k};
        else
            WELLSPECS{i,k} = WELLSPECS_org{2,k};
        end
    end
    WELLSPECS{i,1}=W(i).name;
    WELLSPECS{i,3}=xi(1);
    WELLSPECS{i,4}=xj(1);
    WELLSPECS{i,5}=W(i).refDepth;
end
control.WELSPECS=WELLSPECS;

%% comp dat
%COMPDAT_org=control.COMPDAT;
COMPDAT=cell(0,14);
count=1;
for i=1:numel(W)
    [xi,xj,xk]=ind2sub(cartDims, indexMap(W(i).cells));
    for j=1:numel(W(i).cells)
        for k=1:14
            if(W(i).sign>0)
                COMPDAT{count,k} = COMPDAT_org{1,k};
            else
                COMPDAT{count,k} = COMPDAT_org{2,k};
            end
        end
        COMPDAT{count,1}=W(i).name;
        COMPDAT{count,2}=xi(j);
        COMPDAT{count,3}=xj(j);
        COMPDAT{count,4}=xk(j);
        COMPDAT{count,5}=xk(j);
        if(W(i).status==1)
             COMPDAT{count,6}='OPEN';
        else
            COMPDAT{count,6}='SHUT';
        end
        if ( numel(W(i).r) == numel(W(i).cells) )
            COMPDAT{count,9}=W(i).r(j)*2;
        else
            COMPDAT{count,9}=W(i).r*2;
        end
        if(opt.add_wellindex)
            COMPDAT{count,8}=W(i).WI(j)/((centi*poise*meter^3)/(day*barsa));
        end
        COMPDAT{count,13}=upper(W(i).dir(j));
        count=count+1;
    end
end
control.COMPDAT=COMPDAT;
%%
% wconinje
%WCONINJE_org=deck.SCHEDULE.control(1).WCONINJE;
WCONINJE=cell(0,size(WCONINJE_org,2));
count=1;
for i=1:numel(W)
    if(W(i).sign>0)
        for k=1:size(WCONINJE_org,2)
            WCONINJE{count,k}=WCONINJE_org{1,k};
        end
        WCONINJE{count,1}=W.name;
        if(W(i).compi(:,1)==1)
            itype = 'WATER';
        elseif(W(i).compi(:,2)==1)
            if(RUNSPEC.OIL==1)
                itype = 'OIL';
            else
                itype ='GAS';
            end
        else
            itype = 'GAS';
        end
        WCONINJE{count,1}=W(i).name;
        WCONINJE{count,2}=itype;
        if(W(i).status)
            WCONINJE{count,3}='OPEN';
        else
            WCONINJE{count,3}='SHUT';
        end
        WCONINJE{count,4}=upper(W(i).type);
        switch W(i).type
            case 'bhp'
                WCONINJE{count,7} = W(i).val/barsa;
            case 'resv'
                WCONINJE{count,6} = abs(W(i).val*day);
            case 'rate'
                WCONINJE{count,5} = abs(W(i).val*day);
            otherwise
                error('Not done');
        end 
        if(~isempty(W(i).lims))
            warning('Limits not done');
            %WCONINJE{i,5:7}
        end
        count=count+1;
    end
end

control.WCONINJE=WCONINJE;
%%
% wconprod
% wconinje
%WCONPROD_org=control.WCONPROD;
WCONPROD=cell(0,size(WCONPROD_org,2));
count=1;
for i=1:numel(W)
    if(W(i).sign<0)
        for k=1:size(WCONPROD_org,2)
            WCONPROD{count,k}=WCONPROD_org{1,k};
        end
        WCONPROD{count,1}=W(i).name;
        if(W(i).compi(:,1)==1)
            itype = 'WATER';
        elseif(W(i).compi(:,2)==1)
            if(RUNSPEC.OIL==1)
                itype = 'OIL';
            else
                itype ='GAS';
            end
        else
            itype = 'GAS';
        end
        WCONPROD{count,2}=itype;
        if(W(i).status)
            WCONPROD{count,2}='OPEN';
        else
            WCONPROD{count,2}='SHUT';
        end
        WCONPROD{count,3}=upper(W(i).type);
        switch W(i).type
            case 'bhp'
                WCONPROD{count,9} = W(i).val/barsa;
            case 'resv'
                WCONPROD{count,8} = abs(W(i).val*day);
            case 'orat'
                WCONPROD{count,4} = abs(W(i).val*day);
            case 'grat'
                WCONPROD{count,6} = abs(W(i).val*day);
            case 'wrat'
                WCONPROD{count,5} = abs(W(i).val*day);
            otherwise
                error('Not done');
        end 
                
        if(~isempty(W(i).lims))
            warning('Limits not done');
            %WCONPROD{i,4:13}
        end
        count = count+1;
    end
end
 control.WCONPROD=WCONPROD;
