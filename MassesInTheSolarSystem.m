%%% Inport data for satellite position in orbit

x_axis = [-3.62 -3.46 -3.25 -2.97 -2.60 -2.14 -1.55 -0.85 -0.03 0.78 1.45 1.87 2.09 2.16 2.11 1.99 1.82 1.61 1.37 1.11 0.85 0.58 0.28 0.00 -0.27 -0.56 -0.84 -1.12 -1.38 -1.64 -1.89 -2.14 -2.37 -2.59 -2.80 -2.99 -3.17 -3.33 -3.49 -3.59 -3.69 -3.77 -3.81 -3.83 -3.81 -3.76 -3.65 -3.51];
y_axis = [1.04 0.63 0.20 -0.22 -0.65 -1.03 -1.37 -1.58 -1.59 -1.32 -0.79 -0.11 0.58 1.22 1.82 2.35 2.81 3.22 3.59 3.90 4.16 4.40 4.58 4.74 4.86 4.95 5.01 5.03 5.04 5.00 4.95 4.87 4.77 4.65 4.50 4.33 4.14 3.93 3.69 3.42 3.15 2.85 2.52 2.20 1.83 1.46 1.06 0.65];
time_elapsed = 0:0.25:11.75;

%%% Graph data to show entire orbit

scatter(x_axis, y_axis)
hold on;
scatter(0,0,'black')

%%% Estimate semi-major axis of orbit

[d1,d2] = size(x_axis);
distances = 0:0:d2;

for index = 1:d2
    distances(index) = sqrt(x_axis(index)^2 + y_axis(index)^2);
end

max_dist = max(distances);
min_dist = min(distances);
farthest = [x_axis(find(distances == max_dist)), y_axis(find(distances == max_dist))];
closest = [x_axis(find(distances == min_dist)), y_axis(find(distances == min_dist))];

semi_major_axis = (sqrt((farthest(2)-closest(2))^2 + (farthest(1)-closest(1))^2))/2; % lunar radii

%%% Estimate eccentricity of orbit

% Here, we use the fact that the satellite is basically massless when
% compared with the moon, so the moon is basically the center of mass of
% the system
eccentricity = (max_dist - min_dist)/(max_dist + min_dist);

%%% Estimate mass of moon

% let the period of the orbit be 11.5 hours with 0.14 uncertainty
orbital_period = 11.5 * 60 * 60; % in seconds
mass_of_system = (semi_major_axis*1.7371*10^6)^3*4*pi^2/(6.67408*10^-11*(orbital_period)^2); 