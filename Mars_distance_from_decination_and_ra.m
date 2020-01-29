earth_679day_distance_separation = 0.8825; % in au;
earth_679day_angular_separation = 52.2519;  % in degrees;
earth_sun_distance = 1; % in au;

% days between date 1 and date 2 is exactly 679 days because that is the
% calculated period of mars.

date1_mars_dec = 0;
date1_mars_RA = 0;
date1_sun_dec = 0;
date1_sun_RA = 0;

date2_mars_dec = 0;
date2_mars_RA = 0;
date2_sun_dec = 0;
date2_sun_RA = 0;

earth1_theta = acos(sin(date1_mars_dec)*sin(date1_sun_dec)+cos(date1_mars_dec)*cos(date1_sun_dec)*cos(date1_mars_RA - date1_sun_RA)); % this is the angle of separation between mars and the sun as seen from earth on date1
earth2_theta = acos(sin(date2_mars_dec)*sin(date2_sun_dec)+cos(date2_mars_dec)*cos(date2_sun_dec)*cos(date2_mars_RA - date2_sun_RA)); % this is the angle of separation between mars and the sun as seen from earth on date1

if earth1_theta > 63.875 && earth2_theta > 63.875
     % case 1: one of the thetas is greater than 180
     if earth1_theta > 180 || earth2_theta > 180
         if earth2_theta > 180
             earth2_theta = 360 - earth2_theta;
         end
         if earth1_theta > 180
             earth1_theta = 360 - earth1_theta;
         end
         mars_angle_between_earth1_earth2 = 360 - earth1_theta - earth2_theta - earth_679day_angular_separation; % in degrees
         earth2_dist_to_mars = sin(earth1_theta-63.85)*(0.8825/sin(mars_angle_between_earth1_earth2));
         earth1_dist_to_mars = sin(earth2_theta-63.85)*(0.8825/sin(mars_angle_between_earth1_earth2));
     end
     % case 2: one of the thetas is 180
     if earth1_theta == 180 || earth2_theta == 180
         if earth1_theta == 180
             mars_angle_between_earth1_earth2 = 180 - (earth_679day_angular_separation+earth2_theta);
             earth2_dist_to_mars = sin(earth_679day_angular_separation)*(1/sin(mars_angle_between_earth1_earth2));
             earth1_dist_to_mars = sin(earth2_theta)*(1/sin(mars_angle_between_earth1_earth2));
         end
         if earth2_theta == 180
             mars_angle_between_earth1_earth2 = 180 - (earth_679day_angular_separation+earth1_theta);
             earth1_dist_to_mars = sin(earth_679day_angular_separation)*(1/sin(mars_angle_between_earth1_earth2));
             earth2_dist_to_mars = sin(earth1_theta)*(1/sin(mars_angle_between_earth1_earth2));
         end
     end
     % case 3: both thetas are less than 180
     if earth1_theta < 180 && earth2_theta < 180
          mars_angle_between_earth1_earth2 = 360 - earth1_theta - earth2_theta - earth_679day_angular_separation; % in degrees
          earth2_dist_to_mars = sin(earth1_theta-63.85)*(0.8825/sin(mars_angle_between_earth1_earth2));
          earth1_dist_to_mars = sin(earth2_theta-63.85)*(0.8825/sin(mars_angle_between_earth1_earth2));
     end
end

