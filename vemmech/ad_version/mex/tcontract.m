function varargout = tcontract(varargin)

   filename = 'tcontract.cpp';
   INCLUDE = {};
   OPTS= {'-O'};
   SRC = {filename};
   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();
   
   %@@ DEBUG
   %CXXFLAGS = {'CXXFLAGS=$CXXFLAGS -D_GNU_SOURCE -D_GLIBCXX_DEBUG=1 -g -fPIC -O3 -std=c++11 -ffast-math -march=native -fopenmp'}
   
   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   
   [varargout{1:nargout}] = tcontract(varargin{:});
end

