function SaveVideo(frames,fileName)
        axis tight manual 
        set(gca,'nextplot','replacechildren');
        format = 'Uncompressed AVI';
        v = VideoWriter(fileName,format);
        compressionRatio = 5; 
        frameRate = 5;
        quality = 100;
        v.FrameRate = frameRate;
%         v.CompressionRatio = compressionRatio;
%         v.Quality = quality;
        open(v)
        writeVideo(v,frames);
        close(v)
end