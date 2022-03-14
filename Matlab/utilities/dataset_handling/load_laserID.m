function laserID = load_laserID(path,datasetprefix)
    % from log, retreive laser indices for all scanlines in every scan.
    switch upper(datasetprefix)
        case "A"
            if isfile(strcat(path,filesep,'laserID_log_A.mat'))
                load('laserID_log_A.mat');
            else
                laserID = laserID_logParser_v2([path,'/','laserID_log_A.txt']);
                save([path,'/','laserID_log_A.mat'],'laserID');
            end
        case "P"
            if isfile(strcat(path,filesep,'laserID_log_P.mat'))
                load('laserID_log_P.mat');
            else
                laserID = laserID_logParser_v2([path,'/','laserID_log_P.txt']);
                save([path,'/','laserID_log_P.mat'],'laserID');
            end
        case "SL"
            if isfile(strcat(path,filesep,'laserID_log_sl.mat'))
                load('laserID_log_sl.mat');
            else
                laserID = laserID_logParser_v2([path,'/','laserID_log_sl.txt']);
                save([path,'/','laserID_log_sl.mat'],'laserID');
            end 
    end
end