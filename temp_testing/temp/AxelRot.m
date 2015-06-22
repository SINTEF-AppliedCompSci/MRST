function varargout=AxelRot(varargin)
%Generate roto-translation matrix for the rotation around an arbitrary line in 3D.
%The line need not pass through the origin. Optionally, also, apply this
%transformation to a list of 3D coordinates.
%
%SYNTAX 1:
%
%    M=AxelRot(deg,u,x0)
%
%
%in:
%
%  u, x0: 3D vectors specifying the line in parametric form x(t)=x0+t*u 
%         Default for x0 is [0,0,0] corresponding to pure rotation (no shift).
%         If x0=[] is passed as input, this is also equivalent to passing
%         x0=[0,0,0].
%
%  deg: The counter-clockwise rotation about the line in degrees. 
%       Counter-clockwise is defined using the right hand rule in reference
%       to the direction of u.
%
%
%out:
%
% M: A 4x4 affine transformation matrix representing
%    the roto-translation. Namely, M will have the form
% 
%                 M=[R,t;0 0 0 1] 
%   
%    where R is a 3x3 rotation and t is a 3x1 translation vector. 
% 
%
%
%SYNTAX 2:
%
%       [R,t]=AxelRot(deg,u,x0)
%
% Same as Syntax 1 except that R and t are returned as separate arguments.
% 
%
%
%SYNTAX 3: 
%
% This syntax requires 4 input arguments be specified, 
%
%   [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0)
% 
% where the columns of the 3xN matrix XYZold specify a set of N point
% coordinates in 3D space. The output XYZnew is the transformation of the
% columns of XYZold by the specified rototranslation about the axis. All 
% other input/output arguments are as before.
%
%   by Matt Jacobson
%  
%   Copyright, Xoran Technologies, Inc. 2011


if nargin>3
    
   XYZold=varargin{1};
   varargin(1)=[];
    
   [R,t]=AxelRot(varargin{:});
    
   XYZnew=bsxfun(@plus,R*XYZold,t);
   
   varargout={XYZnew, R,t};
   
   return; 
   
end

    [deg,u]=deal(varargin{1:2});
    
    if nargin>2, x0=varargin{3}; end

    R3x3 = nargin>2 && isequal(x0,'R');

    if nargin<3 || R3x3 || isempty(x0), 
        x0=[0;0;0]; 
    end

    x0=x0(:); u=u(:)/norm(u);

    AxisShift=x0-(x0.'*u).*u;




    Mshift=mkaff(eye(3),-AxisShift);

    Mroto=mkaff(R3d(deg,u));

    M=inv(Mshift)*Mroto*Mshift;

    varargout(1:2)={M,[]};
    
    if R3x3 || nargout>1 
      varargout{1}=M(1:3,1:3);
    end
    
    if nargout>1,
      varargout{2}=M(1:3,4);  
    end

    
    
function R=R3d(deg,u)
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%
%R=R3d(deg,axis)
%
% deg: The counter-clockwise rotation about the axis in degrees.
% axis: A 3-vector specifying the axis direction. Must be non-zero

    R=eye(3);
    u=u(:)/norm(u);
    x=deg; %abbreviation

    for ii=1:3

        v=R(:,ii);

        R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
          %Rodrigues' formula

    end



function M=mkaff(varargin)
% M=mkaff(R,t)
% M=mkaff(R)
% M=mkaff(t)
%
%Makes an affine transformation matrix, either in 2D or 3D. 
%For 3D transformations, this has the form
%
% M=[R,t;[0 0 0 1]]
%
%where R is a square matrix and t is a translation vector (column or row)
%
%When multiplied with vectors [x;y;z;1] it gives [x';y';z;1] which accomplishes the
%the corresponding affine transformation
%
% [x';y';z']=R*[x;y;z]+t
%



    if nargin==1

       switch numel(varargin{1}) 

           case {4,9} %Only rotation provided, 2D or 3D

             R=varargin{1}; 
             nn=size(R,1);
             t=zeros(nn,1);

           case {2,3}

             t=varargin{1};
             nn=length(t);
             R=eye(nn); 

       end
    else

        [R,t]=deal(varargin{1:2});
        nn=size(R,1);
    end

    t=t(:); 

    M=eye(nn+1);

    M(1:end-1,1:end-1)=R;
    M(1:end-1,end)=t(:); 

