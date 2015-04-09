%%



mrstVerbose off
t = tic;
upscaler = TwoPhaseUpscaler(G, rock, fluid);
upscaler.partition = p;
%upscaler.blocks = 2;
upscaler.dims = 1:3;
upscaler.nvalues = 19;
upscaler.method = 'capillary';
dataTP = upscaler.upscale();
toc(t)


