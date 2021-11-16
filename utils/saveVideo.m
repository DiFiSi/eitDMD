function saveVideo(vid, filename, fr, cmap)
    vid = double(vid);
    
    if ~exist('fr','var')
       fr = 15;
    end
    
    if ~exist('cmap','var')
       cmap = 'jet';
    end
    
    v = VideoWriter(filename);
    v.FrameRate = fr;
    
    open(v);
    
    figure;
    for s = 1:size(vid,3) 
        imagesc(vid(:,:,s)); colormap(cmap);
        frame = getframe(gcf);
        
        writeVideo(v,frame);
    end
    
    close all;
    close(v);
end