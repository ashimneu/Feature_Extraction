function clearfig(fignum)
    % If exists, clear the figure whose number is fignum.
    
    % get handle of figures for the specified figure number
    fighandle = findobj('type','figure','Number',fignum); 

    % clear the figure if it exists.
    if ~(isempty(fighandle))
        clf(fignum,'reset');        
    end    
end