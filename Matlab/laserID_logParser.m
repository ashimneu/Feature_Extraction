function Ids =  laserID_logParser(path)
Ids= struct;
    fid = fopen(path);
    
    count = 1;
while ~feof(fid)
    tline = fgetl(fid);
    ID_text = strsplit(tline,',');
    Ids(count).ID = str2num( cell2mat(ID_text));
    count = count + 1;
end

    

end