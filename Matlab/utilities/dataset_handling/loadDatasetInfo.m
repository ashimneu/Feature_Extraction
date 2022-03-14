function [path,scanlist,datasetprefix] = loadDatasetInfo(dataset_number)
    if dataset_number == 1
            % Dataset 1 - ESS dataset
            path = 'data/data_8LASERS'; 
            scanlist = dir([path,'/wrong_est/*.pcd']); % path of the source folder of scans
            datasetprefix = "P";
    elseif dataset_number == 2
            % Dataset 2 - ESS dataset with range & theta, use prefix - P
            path = 'data/data_8LASERS_ess_20211222120852'; 
            scanlist = dir([path,'/wrong_est/*.pcd']); % path of the source folder of scans
            datasetprefix = "P";
    elseif dataset_number == 3
            % Dataset 3- Dataset of container yard
            path = 'data/container_dataset_7_20ft_17pan'; % 
            scanlist = dir([path,'/scanlines7_20ft_pan17/*.pcd']); % path of the source folder of scans
            datasetprefix = "sl";
    else
        error('The dataset you chose doesn''t exist.')

    end
end

