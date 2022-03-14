function clearaxis(fignum)
    % If exists, clear axis of the figure whose fig. number = fignum.
    
    % get handle of figures for the specified figure number
    fighandle = findobj('type','figure','Number',fignum); 

    % clear the figure if it exists.
    if ~(isempty(fighandle))
        cla;        
    end    
end