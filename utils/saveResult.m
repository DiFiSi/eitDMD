function saveResult(dataFile, data, fig)
    [parentFolder,~,~] = fileparts(dataFile);
    
    if ~exist(parentFolder, 'dir')
        mkdir(parentFolder);
    end
    
    save(dataFile, 'data');
    savefig(fig, dataFile);
    % insert save .pdf
end