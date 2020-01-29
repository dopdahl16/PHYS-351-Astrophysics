% Create Salpeter distribution of 400 stars

solar_mass = 1.988 * 10^30; % kg
rand_dist = []; % array of 400 random numbers 0-1
tot_num_stars = 400;
for i = 1:tot_num_stars
   rand_dist(i) = rand;
end
solar_masses = 0.5:0.1:30.1; % In terms of solar masses - 2 means 2 * solar_mass
rand_dist = sort(rand_dist);
D = 1/((30*solar_mass)^-1.35 - (0.5*solar_mass)^-1.35);
for index = 1:length(solar_masses)-1
    high = solar_masses(1,index+1);
    low = solar_masses(1,index);
    bin_weight = D*((high*solar_mass)^-1.35 - (low*solar_mass)^-1.35);
    solar_masses(2,index) = bin_weight;
end
for index = 1:length(solar_masses)
    upper_limit = 0;
    for h = 1:index
        upper_limit = upper_limit + solar_masses(2,h);
    end
    count = 0;
    for i = 1:length(rand_dist)
        if rand_dist(i) <= upper_limit
            count = count + 1;
        end
    end
    solar_masses(3,index) = count;
end
solar_masses(:,297) = [];
previous_val = 0;
for j = 1:length(solar_masses)
     g = solar_masses(3,j);
     solar_masses(3,j) = solar_masses(3,j) - previous_val;
     previous_val = g;
end
% Assign 1/3 of stars a binary partner
binary_solar_masses = solar_masses(1:2,:);
binary_solar_masses(3,1) = 0;
limit = 133;
for bi_index = 1:limit
    bi_star = round(rand()*399 + 1);
    count = 0;
    for index = 1:length(solar_masses)
        for h = 1:solar_masses(3,index)
            count = count + 1;
            if count == bi_star
                if solar_masses(3,index) == binary_solar_masses(3,index)
                    limit = limit + 1;
                else
                    binary_solar_masses(3,index) = binary_solar_masses(3,index) + 1;
                end
            end
        end
    end
end

for index = 1:length(solar_masses)
    solar_masses(3,index) = solar_masses(3,index) - binary_solar_masses(3,index);
end
for index = 1:length(solar_masses)
    if solar_masses(3,index) == -1
        disp(index)
    end
end   
tot = 0;
bi_tot = 0;
for index = 1:length(solar_masses)
    tot = solar_masses(3,index) + tot
    bi_tot = binary_solar_masses(3,index) + bi_tot
end  
% Create ZAMS HR diagram
scatter(solar_masses(1,:),solar_masses(3,:))