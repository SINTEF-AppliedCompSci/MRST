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
  
3. Line 90 of DFN.m (pts = rand(n,2);) needs to be changed to:

         if isempty(opt.pos)
         
                  pts = [Scale(rand(n,1),opt.bbx(1),opt.bbx(3)),Scale(rand(n,1),...   % location bbox(0,1)
                  
                  opt.bbx(2),opt.bbx(4))];
                  
         else
         
                  pts = opt.pos;
                  
         end
   
4. Lines 167 - 172 of DFN.m needs to be changed to:

         n = opt.n;
         
         if isempty(opt.pos)                                                     % (!) added box limits
        
                  pts = [Scale(rand(n,1),opt.bbx(1),opt.bbx(4)),...                   % bounded U for each axis
            
                  Scale(rand(n,1),opt.bbx(2),opt.bbx(5)),...
                
                  Scale(rand(n,1),opt.bbx(3),opt.bbx(6))];
                
         else
        
                  pts = opt.pos;                                                      % any other desired points
            
         end
        
         T = [1,0,0,pts(i,1);0,1,0,pts(i,2);0,0,1,pts(i,3);0,0,0,1]*...                % transformers
        
                  [opt.scale,0,0,0;0,opt.scale,0,0;0,0,opt.scale,0;0,0,0,1]*...       % scaling
            
                  [cos(dds(i)),-sin(dds(i)),0,0;sin(dds(i)),cos(dds(i)),0,0;...       % rotations
            
                  0,0,1,0;0,0,0,1]*[cos(dps(i)),0,sin(dps(i)),0;0,1,0,0;...
            
                  -sin(dps(i)),0,cos(dps(i)),0;0,0,0,1];
                  
                           
