triton_vectors = [[52227239.9031555,34074967.1972109,35741068.0887101,43423629.8251168,42831491.2365865,47987518.2603380,60773831.4202242,56752221.3000091;0.442883820350673,-1.45646458147664,0.371001712840482,-0.399308222933237,1.39964772137871,-0.485226992567337,0.681770703400577,-0.448723344010721;22.9200000000000,46.9700000000000,94.8800000000000,142.650000000000,166.570000000000,190.500000000000,214.150000000000,238.350000000000]];
% for i = 1:8
%     
% TITLE = 'Select a .jpg image to calibrate distances:';
% [cal_filename,cal_pathname] = uigetfile('*.*',TITLE,'C:\','MultiSelect','off'); 
% 
% cal_full_filename = [cal_pathname cal_filename];
% current_image = imread(cal_full_filename);
% inverted_current_image = 255-current_image;
% 
% image(inverted_current_image,'CDataMapping','scaled');
% 
%     colormap gray
%     axis image
%     title('Calibration image now showing')
%     drawnow;
% 
%     %calibration based off of neptunes size
% 
%     [x1_cal, y1_cal, ~] = ginput(1);
% 
%     [x2_cal, y2_cal, ~] = ginput(1);
% 
%     x_scale = 49.24e6 / abs(x2_cal - x1_cal);
% 
%     [x3_cal, y3_cal, ~] = ginput(1);
% 
%     [x4_cal, y4_cal, ~] = ginput(1);
% 
%     y_scale = 49.24e6 / abs(y3_cal - y4_cal);
% 
%     [x1,y1,~] = ginput(1);
% 
%     [x2,y2,~] = ginput(1);
%     
%     x_dist = (x1 - x2) * x_scale;
%     y_dist = (y1 - y2) * y_scale;
%     
%     magnitude = sqrt(x_dist^2 + y_dist^2);
%     angle = atan(y_dist/x_dist);
%     
%     triton_vectors(1,i) = magnitude;
%     triton_vectors(2,i) = angle;
%     
% end
% triton_vectors

% add the times to the matrix triton_vectors in hours
% 0 was midnight on Aug 28, 2002

triton_vectors(3,1) = 22.92;
triton_vectors(3,2) = 46.97;
triton_vectors(3,3) = 94.88;
triton_vectors(3,4) = 142.65;
triton_vectors(3,5) = 166.57; 
triton_vectors(3,6) = 190.5;
triton_vectors(3,7) = 214.15;
triton_vectors(3,8) = 238.35;

% to determine the semimajor axis, we average the magnitude vectors

semi_major_axis = (triton_vectors(1,1)+triton_vectors(1,2)+triton_vectors(1,3)+triton_vectors(1,4)+triton_vectors(1,5)+triton_vectors(1,6)+triton_vectors(1,7)+triton_vectors(1,8))/8
t = triton_vectors';
unc_semi_major_axis = std(t)

% to determine the period, I selected two angles that were the closest to
% eachother and found the difference in time between the two

angular_velocity = 1:7;

for i = 1:7
    time_elapssed = (triton_vectors(3, i+1) - triton_vectors(3, i)) * 60 * 60;
    angle_covered = triton_vectors(2, i+1) - triton_vectors(2,i);
    if angle_covered < 0
        angle_covered = angle_covered + 2*pi;
    end
    angular_velocity(i) = angle_covered / time_elapssed;
end

avg_angular_speed = mean(angular_velocity)
unc_avg_angular_speed = std(angular_velocity)

period = (avg_angular_speed / (2*pi))^(-1)

mass_neptune = (semi_major_axis^3*4*pi^2)/(6.67408*10^-11*period^2)
unc_mass_neptune = sqrt((unc_semi_major_axis(1))^2 + (unc_avg_angular_speed)^2)
