function sm_xyz = doSmoothing_Plot(pc_xyz,sm_window)
 % Performs smoothing of data using simple moving average on the data
 % average is computed in a neighboorhood of a point using a sliding window
 
 sm_xyz = smoothdata(pc_xyz,'movmean',sm_window); % smoothened points
 sz_marker = 6; % marker size
 figure(7); clf; hold on; grid on;
 scatter3(sm_xyz(:,1),sm_xyz(:,2),sm_xyz(:,3),sz_marker,'b');
 xlabel('x');
 ylabel('y');
 zlabel('z');

end

