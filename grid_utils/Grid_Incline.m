function [Grid] = Grid_Incline(Lx,Ly,Nx,Ny,iangle)
%Generate inclined grid
%   iangle - internal inclined angle
angle_limit=atan(Ly/Lx);
if(iangle<angle_limit||iangle>(pi-angle_limit)||iangle==pi/2)
    error('Invalid incline angle!')
else
    Grid=cartGrid([Nx Ny],[Lx Ly]);
    
    Grid.cells.nodes=zeros(Grid.cells.num,4);
    for i=1:Grid.cells.num
        myfaces=Grid.cells.faces(Grid.cells.facePos(i):...
            Grid.cells.facePos(i+1)-1);
        mynodes=cell(numel(myfaces),1);
        for j=1:numel(myfaces)
            myface=myfaces(j);
            mynodes(j)={Grid.faces.nodes(Grid.faces.nodePos(myface):...
                Grid.faces.nodePos(myface+1)-1)};
        end
        mynodes=cell2mat(mynodes);
        mynodes=unique(mynodes);
        x=Grid.nodes.coords(mynodes,1);
        y=Grid.nodes.coords(mynodes,2);
        [~,ind]=min(x.^2+y.^2);
        Grid.cells.nodes(i,1)=mynodes(ind);x(ind)=[];y(ind)=[];mynodes(ind)=[];
        [~,ind]=max(x.^2+y.^2);
        Grid.cells.nodes(i,3)=mynodes(ind);x(ind)=[];y(ind)=[];mynodes(ind)=[];
        [~,ind]=max(x);Grid.cells.nodes(i,2)=mynodes(ind);
        [~,ind]=max(y);Grid.cells.nodes(i,4)=mynodes(ind);
    end
    
    coords=Grid.nodes.coords;
    coords(:,1)=coords(:,1)-0.5*Lx;
    coords(:,2)=coords(:,2)-0.5*Ly;
    
    for i=1:size(coords,1)
        xi=coords(i,1);
        if(xi~=-0.5*Lx&&xi~=0.5*Lx)
            sequence=-Nx/2+1:Nx/2-1;
            for j=1:length(sequence)
                xcoor=Lx/Nx*sequence(j);
                if(abs(xi-xcoor)<1e-10)
                    break;
                end
            end
            if(xcoor<0)
                slope=tan(iangle);
                d1=abs(-0.5*Lx+0.5*Ly/slope);
                d2=abs(-0.5*Lx-0.5*Ly/slope);
                d1=d1/0.5/Nx;
                d2=d2/0.5/Nx;
                x1=-0.5*Lx+j*d1; y1=-0.5*Ly;
                x2=-0.5*Lx+j*d2; y2=0.5*Ly;
                temp=(y1-y2)/(x1-x2);
                coords(i,1)=x1+(coords(i,2)-y1)/temp;
            elseif(xcoor>0)
                slope=tan(iangle);
                d1=abs(-0.5*Ly/slope-0.5*Lx);
                d2=abs(0.5*Ly/slope-0.5*Lx);
                d1=d1/0.5/Nx;
                d2=d2/0.5/Nx;
                x1=-0.5*Ly/slope+(j-0.5*Nx)*d1; y1=-0.5*Ly;
                x2=0.5*Ly/slope+(j-0.5*Nx)*d2;y2=0.5*Ly;
                temp=(y1-y2)/(x1-x2);
                coords(i,1)=x1+(coords(i,2)-y1)/temp;
            else
                coords(i,1)=coords(i,2)/tan(iangle);
            end
        end
    end
    Grid.nodes.coords=coords;
    Grid=computeGeometry(Grid);
end
end
