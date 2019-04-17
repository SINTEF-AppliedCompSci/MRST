function [G] = importGrid(filename,opt)
% import triangle or tetrahedral grid
% opt - 'tri'
%     - 'tetBench'

fid=fopen(filename);
switch opt
    case 'tri'  %---------------------------------------------------------
        while ~feof(fid)
            myline=fgetl(fid);
            if(strfind(myline,'vertices'))
                myline=fgetl(fid);
                np=textscan(myline,'%d');
                np=np{1};
                p=zeros(np,2);
                break;
            end
        end
        for i=1:np
            myline=fgetl(fid);
            data=textscan(myline,'%f %f');
            p(i,1)=data{1};p(i,2)=data{2};
        end
        while ~feof(fid)
            myline=fgetl(fid);
            if(strfind(myline,'triangles'))
                myline=fgetl(fid);
                nt=textscan(myline,'%d');
                nt=nt{1};
                t=zeros(nt,3);
                break;
            end
        end
        for i=1:nt
            myline=fgetl(fid);
            data=textscan(myline,'%d %d %d');
            t(i,1)=data{1};t(i,2)=data{2};t(i,3)=data{3};
        end
        fclose(fid);
        G=triangleGrid(p,t);
    case 'tetBench' %------------------------------------------------------
        while ~feof(fid)
            myline=fgetl(fid);
            if(strfind(myline,'Vertices'))
                np=textscan(myline,'%s %f');
                np=np{2};
                p=zeros(np,3);
                break
            end
        end
        i=1;
        while ~feof(fid)
            myline=fgetl(fid);
            if(isempty(strfind(myline,'Volumes->faces')))
                data=textscan(myline,'%f %f %f');
                p(i,1)=data{1};p(i,2)=data{2};p(i,3)=data{3};i=i+1;
            else
                break;
            end
        end       
        while ~feof(fid)
            myline=fgetl(fid);
            if(strfind(myline,'Volumes->Verticess'))
                nt=textscan(myline,'%s %f');
                nt=nt{2};
                t=zeros(nt,4);
                break;
            end
        end       
        i=1;
        while ~feof(fid)
            myline=fgetl(fid);
            if(isempty(strfind(myline,'Faces->Edgess')))
                data=textscan(myline,'%f %f %f %f %f');
                t(i,1)=data{2};t(i,2)=data{3};t(i,3)=data{4};t(i,4)=data{5};i=i+1;
            else
                break;
            end
        end  
        fclose(fid);
        G=tetrahedralGrid(p,t);
end
end

