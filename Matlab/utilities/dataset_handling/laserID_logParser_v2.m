function Ids =  laserID_logParser_v2(path)
    Ids= struct;
    fid = fopen(path);    
    count = 1;
    while ~feof(fid)
        tline = fgetl(fid);
        line_char    = strsplit(tline,','); % each line is read in as string. % line entries: scan index, laser ID, range, theta, laser ID, range, theta, ... 
        line_string  = convertCharsToStrings(line_char);
        line_double = str2double(line_string); % convert text line to array of numeric values
        Ids(count).scan_idx = line_double(1);       % keep scan index only
        Ids(count).ID       = line_double(2:3:end); % keep laser IDs only
        Ids(count).range    = line_double(3:3:end); % keep range only
        Ids(count).theta    = line_double(4:3:end); % keep theta only
        count = count + 1;
    end
end