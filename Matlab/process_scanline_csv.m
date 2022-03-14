%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ReadMe  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script read in the saved scans from the C++ Program and visualize them in Matlab as an 
% basic frame to start works on extracting features from single scans. 
% Please modify the parameters as wished before running the code.
% Parameters: 
%     idx_N: The script currently uses one lasers if there are multiple lasers in the scan. 
% This parameter figure out which to use.Be careful with the index. e.g. The laser index of M8 is
% from 0-7 for the 8 lasers. While in Matlab, the first index is 1, not 0.
%     FLAG_show_Accum: Debugging flag, if you want to see the accumulation process through the 
% scans, turn this on. But be awared that there is no pre-allocation of memory so this will run
% slower and slower when the point cloud is accumulated.
%     FLAG_csv_or_mat: Data source flag. If you have the .mat data file, you can set this as 1 
% to read data from it. Or if you have the .csv data files for the scans. You can choose the 
% folder containing the .csv files to parse them and run.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
idx_N = 6;                      % Laser ID chosen to analyse
FLAG_show_Accum = 0;            % flag to turn on/off displaying accumulated pointcloud
FLAG_csv_or_mat = 1;            % Choose data source, 0: from CSV, 1: from Mat

%% Start
if FLAG_csv_or_mat
    load("BigShipScans.mat");
else
    selpath = uigetdir;             % Choose the folder saving *.csv form scans
    list = dir(selpath+"\*.csv");   % path of the source folder of scans
end
laser_ang_v = [
    -0.318505;
    -0.2692;
    -0.218009;
    -0.165195;
    -0.111003;
    -0.0557982;
    0;
    0.0557982];             % Manufactor manual shows the vertical angle of each laser in radian
PC_P = [];
for i= 1:length(list) % read through each file in the folder of scans
    if FLAG_csv_or_mat
        scan_table = scans(i).scan_table;                                                           % Read one scan from Mat file
    else
        scan_table = readtable(list(i).folder+"\"+list(i).name,"VariableNamingRule","preserve");    % Read one scan from csv file
    end
    scan = cell2mat(table2cell(scan_table)); % convert the scan to a matrix form. Refer to the table for definition of the fields

    %% Extract point from laser-N
    idx = find(scan(:,2)==idx_N);
    scanline = scan(idx,:);

    %% Convert from (Range, horizontal angle, vertical angle) to L-frame (x,y,z) coordinates
    ver_ang = laser_ang_v(idx_N+1);             % Get vertical angle of this laser, plus 1 due to change of indexing
    scan_L = zeros(size(scanline,1),3);         % allocate memory
    for j=1:size(scan_L,1)
        scan_L(j,1) = scanline(j,3) * cos(ver_ang) * cos(scanline(j,4)); % x = R * cos(θ) * cos(ψ)
        scan_L(j,2) = scanline(j,3) * cos(ver_ang) * sin(scanline(j,4)); % y = R * cos(θ) * sin(ψ)
        scan_L(j,3) = scanline(j,3) * sin(ver_ang);                  % z = R * sin(θ)
    end

    %% Getting saved P-frame coordinates
    scan_P = scan(:,7:9);

    %% plotting
    figure(1)
    subplot(1,2,1)
    scatter3(scan_L(:,1),scan_L(:,2),scan_L(:,3),'.')
    title('L-frame scan calculated from Raw data')
    xlabel('x, m');ylabel('y, m');zlabel('z, m')
    subplot(1,2,2)
    scatter3(scan_P(idx,1),scan_P(idx,2),scan_P(idx,3),'.')
    title('P-frame saved from SHU for comparison')
    xlabel('x, m');ylabel('y, m');zlabel('z, m')

    if FLAG_show_Accum
        % Accumulated points in the P-frame and show, Warning: may delay your PC
        figure(2)
        PC_P = [PC_P;scan_P];
        pcshow(PC_P)
    end

    drawnow
    pause(0.1)
end

