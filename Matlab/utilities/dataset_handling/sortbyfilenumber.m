function list_output = sortbyfilenumber(scanlist,prefix)
    % .pcd list are randomly imported into the array list.
    % Sort in ascending order based on filename.
    % prefix - given non-integer character prefix of .pcd files  

    %     % extract subdirectory name to get prefix
    %     path1_raw = list(1).folder;
    %     path1_split = split(path1_raw,filesep);
    %     subdir = string(path1_split{end}); % dataset subdirectory name
    
    % get filenames & extract integer characters only for sorting
    names_raw = [scanlist.name]; % get filenames
    names_split = split(names_raw,[".pcd",strcat(prefix,"_")]); % remove non-integer characters
    names_str = string(names_split); % convert array of cells to string array
    idx   = names_str == "";   % find empty entries
    names_str(idx) = [];             % remove empty entries
    names = str2double(names_str);   % convert file number to double
    
    % sort in ascending order
    [~, asc_idx] = sort(names,'ascend'); % get index for sorting
    list_output = scanlist(asc_idx);
end