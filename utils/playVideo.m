function h = playVideo(vid, fr, colormap, rangeMin, rangeMax)
    % vid - a 3D matrix (2D over time)
    % fr - framerate
    % colormap - colormap string, e.g. 'autumn' or 'grayscale' or 'jet'
    % rangeMin, rangeMax - lower and higher limits of pixel values in
    % matrix
    
    % make sure data is double
    vid = double(vid);
    
    if ~exist('fr','var')
       fr = 30;
    end
    
    if ~exist('colormap','var')
       colormap = 'jet';
    end
    
    if ~exist('rangeMin','var')
       rangeMin = min(vid(:));
    end
    
    if ~exist('rangeMax','var')
       rangeMax = max(vid(:));
    end
    
    h = implay(vid,fr);
    
    h.Visual.ColorMap.UserRange = 1;
    h.Visual.ColorMap.UserRangeMin = rangeMin;
    h.Visual.ColorMap.UserRangeMax = rangeMax;
    h.Visual.ColorMap.MapExpression = colormap;
    
%     set(0,'showHiddenHandles','on')
%     handle = gcf ;  
%     handle.findobj % to view all the linked objects with the vision.VideoPlayer
%     ftw = handle.findobj ('TooltipString', 'Maintain fit to window');   % this will search the object in the figure which has the respective 'TooltipString' parameter.
%     ftw.ClickedCallback()  % execute the callback linked with this object
end