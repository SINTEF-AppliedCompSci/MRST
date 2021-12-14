The shale module was developed to work with or without the ADFNE package. However, to use ADFNE to generate stochastic fractures in a box with dimensions other than one in all directions, a couple of changes are required.
The changes needed in the DFN.m code are summarized below:

1. Lines 47-49 of DFN.m need to be changed to: 
opt = Option(varargin,'n',100,'minl',0.05,'mu',0.3,'maxl',0.6,...               % default arguments

         'bbx',[0,0,1,1],'polyDomain',{},'dim',2,'asep',0,'dsep',0,'mit',100,'scale',1,...
         
         'shape','c','q',24,'dip',45,'ddip',-1e-7,'dir',0,'ddir',-1e-7,'pos',[]);  % domain added
      
2. Line 64 of DFN.m (cnt = rand(1,2);) needs to be changed to:
cnt = [Scale([0 rand 1],opt.bbx(1),opt.bbx(3)),...                  % random points

         Scale([0 rand 1],opt.bbx(2),opt.bbx(4))];
         
         cnt = [cnt(2) cnt(5)];
  
3. 
