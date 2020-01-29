% Determine luminosicty, temp, mass relationship
mass_temp = [37 23 17.5 7.6 5.9 3.8 2.9 2.0 1.6 1.4 1.05, 1.00 0.79 0.67 0.51 0.4 0.21; 39500 35800 30000 18800 15200 11400 9800 8190 7300 6650 5940 5777 5150 4410 3840 3520 3170]; % in K
mass_lum = [37 23 17.5 7.6 5.9 3.8 2.9 2.0 1.6 1.4 1.05, 1.00 0.79 0.67 0.51 0.4 0.21; 324000 147000 32500 1580 480 96.7 39.4 12.3 5.21 2.56 1.25 1 0.552 0.216 0.077 0.032 0.0076]; % L/L(sol)
mass_temp_lum = [37 23 17.5 7.6 5.9 3.8 2.9 2.0 1.6 1.4 1.05, 1.00 0.79 0.67 0.51 0.4 0.21; 39500 35800 30000 18800 15200 11400 9800 8190 7300 6650 5940 5777 5150 4410 3840 3520 3170; 324000 147000 32500 1580 480 96.7 39.4 12.3 5.21 2.56 1.25 1 0.552 0.216 0.077 0.032 0.0076];
% figure
% scatter(mass_temp(1,:),mass_temp(2,:))
% hold on
f = fit(mass_temp(1,:)',mass_temp(2,:)','linearinterp');
% plot(f)
% title("Mass / Temp");'linearinterp'

% figure
% scatter(mass_lum(1,:),mass_lum(2,:))
% hold on
myfit = fit(mass_lum(1,:)',mass_lum(2,:)','linearinterp');
% plot(myfit)
% title("Mass / Lum");

% figure
% scatter(mass_temp_lum(2,:)', mass_temp_lum(3,:)')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca,'Xdir','reverse')



% Create Salpeter distribution of 400 stars

solar_mass = 1.988 * 10^30; % kg
rand_dist = []; % array of 400 random numbers 0-1
tot_num_stars = 900;
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

% Assign binary partners
star_index = 1:tot_num_stars;
paired_stars = zeros(5,round(tot_num_stars/3));
for index = 1:(round(tot_num_stars)/3)
    bi_star1 = round(rand()*(length(star_index)-1) + 1);
    paired_stars(1,index) = star_index(bi_star1);
    star_index(bi_star1) = [];
    bi_star2 = round(rand()*(length(star_index)-1) + 1);
    paired_stars(2,index) = star_index(bi_star2);
    star_index(bi_star2) = [];
end

% The numbers in paired_stars represent the star's number. Not the mass.
% For example, star 51 is the 51st star that can be found by determining
% where the 51st star is in solar_masses via a running total.
for star_pair = 1:length(paired_stars)
    count = 1;
    for mass_bin = 1:length(solar_masses)
        for star = 1:solar_masses(3,mass_bin)
            if count == paired_stars(1,star_pair)
                paired_stars(4,star_pair) = solar_masses(1,mass_bin);
            end
            count = count + 1;
        end
    end
end
for star_pair = 1:length(paired_stars)
    count = 1;
    for mass_bin = 1:length(solar_masses)
        for star = 1:solar_masses(3,mass_bin)
            if count == paired_stars(2,star_pair)
                paired_stars(5,star_pair) = solar_masses(1,mass_bin);
            end
            count = count + 1;
        end
    end
end
for solo_star = 1:length(star_index)
    count = 1;
    for mass_bin = 1:length(solar_masses)
        for star = 1:solar_masses(3,mass_bin)
            if count == star_index(1,solo_star)
                star_index(3,solo_star) = solar_masses(1,mass_bin);
            end
            count = count + 1;
        end
    end
end
solar_bins_mass_temp_lum = [];
solar_bins_mass_temp_lum(1,:) = solar_masses(1,:);


for index = 1:length(solar_bins_mass_temp_lum)
    mass = solar_bins_mass_temp_lum(1,index);
    % Find the temp based on the mass
    for i = 1:length(mass_temp)-1
        if mass > mass_temp(1,i+1) && mass < mass_temp(1,i)
            slope = (mass_temp(2,i+1)-mass_temp(2,i))/(mass_temp(1,i+1)-mass_temp(1,i));
            intercept = mass_temp(2,i)-(slope*mass_temp(1,i));
            solar_bins_mass_temp_lum(2,index) = slope*solar_bins_mass_temp_lum(1,index)+intercept;
        elseif mass == mass_temp(1,i)
            solar_bins_mass_temp_lum(2,index) = mass_temp(2,i);
        end
    end
    % Find the luminosity based on the mass
    for i = 1:length(mass_lum)-1
        if mass > mass_lum(1,i+1) && mass < mass_lum(1,i)
            slope = (mass_lum(2,i+1)-mass_lum(2,i))/(mass_lum(1,i+1)-mass_lum(1,i));
            intercept = mass_lum(2,i)-(slope*mass_lum(1,i));
            solar_bins_mass_temp_lum(3,index) = slope*solar_bins_mass_temp_lum(1,index)+intercept;
        elseif mass == mass_lum(1,i)
            solar_bins_mass_temp_lum(3,index) = mass_lum(2,i);
        end
    end
end

% for t = 1:length(paired_stars)
%     if paired_stars(4,t) <= paired_stars(5,t)
%         ratio = paired_stars(4,t)/paired_stars(5,t);
%     else
%         ratio = paired_stars(5,t)/paired_stars(4,t);
%     end
%     paired_stars(6,t) = ratio;
% end
% figure
% histogram(paired_stars(6,:))









kb = 8.61734e-5;
c = 3*10^17;
hc = 1240;
e = 2.71828;
b_sigma = 47;%e-9;
b_mu = 445;%e-9;
v_sigma = 44;%e-9;
v_mu = 551;%e-9;
T = 0;

% syms Planck(x)
% Planck = @(x) (2*hc*c/x^5)*(1/(e^(hc/(x*kb*T))-1));
% figure
% fplot(Planck, [0 2000]);

% syms B(x)
% B = @(x) (1/(b_sigma*sqrt(2*pi)))*(exp(-0.5*((x-b_mu)/(b_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T))-1));
% figure
% fplot(B, [0 2000]);
% syms V(x)
% V = @(x) (1/(v_sigma*sqrt(2*pi)))*(exp(-0.5*((x-v_mu)/(v_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T))-1));
% b_integral = integral(B, 0, Inf);
% v_integral = integral(V, 0, Inf);
% 
% bv_ratio = b_integral / v_integral;






temp_bv = [42000 39500 37500 35800 30000 25400 20900 18800 15200 13700 12500 11400 10500 9800 9400 9020 8190 7600 7300 7050 6650 6250 5940 5790 5777 5310 5150 4990 4690 4540 4410 4150 3840 3660 3520 3400 3290 3170 3030 2860; -.33 -.33 -.32 -.32 -.3 -.26 -.24 -.2 -.17 -.15 -.13 -.11 -.07 -.02 .01 .05 .15 .25 .3 .35 .44 .52 .58 .63 .65 .74 .81 .86 .96 1.05 1.15 1.33 1.4 1.46 1.49 1.51 1.54 1.64 1.73 1.8];
% figure
% scatter(temp_bv(1,:), temp_bv(2,:))
% hold on
bvfittype = fittype('a+b*x^n','dependent', {'y'}, 'independent',{'x'},'coefficients', {'a','b','n'});
bvfit = fit(temp_bv(1,:)', temp_bv(2,:)',bvfittype,'StartPoint',[-0.0001,2.0743,-0.0001]);
bvcoeff = coeffvalues(bvfit);
% plot(bvfit)


final_data_single_stars = ["mass";"b-v";"luminosity"];
for solo_star = 1:length(star_index)
    for solar_mass = 1:length(solar_bins_mass_temp_lum)
        T = solar_bins_mass_temp_lum(2,solar_mass);
        if solar_bins_mass_temp_lum(1,solar_mass) == star_index(3, solo_star)
            final_data_single_stars(1,solo_star+1) = star_index(3, solo_star);
            final_data_single_stars(2,solo_star+1) = bvcoeff(1)+bvcoeff(2)*T^bvcoeff(3);
            final_data_single_stars(3,solo_star+1) = solar_bins_mass_temp_lum(3,solar_mass);
        end
    end
end


final_data_binary_stars = ["mass1";"mass2";"b-v1";"b-v2";"b-v_combined";"luminosity1";"luminosity2";"luminosity_combined"];
for star_pair = 1:length(paired_stars)
    final_data_binary_stars(1,star_pair+1) = paired_stars(4, star_pair);
    final_data_binary_stars(2,star_pair+1) = paired_stars(5, star_pair);
    if str2double(final_data_binary_stars(2,star_pair+1)) == 0 || str2double(final_data_binary_stars(1,star_pair+1)) == 0
        disp("WTF")
    end
    for solar_mass = 1:length(solar_bins_mass_temp_lum)
        if solar_bins_mass_temp_lum(1,solar_mass) == paired_stars(4, star_pair)
            T1 = solar_bins_mass_temp_lum(2,solar_mass);
            Bx1 = @(x) (1/(b_sigma*sqrt(2*pi)))*(exp(-0.5*((x-b_mu)/(b_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T1))-1));
            Vx1 = @(x) (1/(v_sigma*sqrt(2*pi)))*(exp(-0.5*((x-v_mu)/(v_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T1))-1));
            b_integral1 = integral(Bx1, 300, 1000);
            v_integral1 = integral(Vx1, 300, 1000);
            bv_ratio1 = b_integral1 / v_integral1;
            final_data_binary_stars(3,star_pair+1) = bvcoeff(1)+bvcoeff(2)*solar_bins_mass_temp_lum(2,solar_mass)^bvcoeff(3);
            final_data_binary_stars(6,star_pair+1) = solar_bins_mass_temp_lum(3,solar_mass);
        end
        if solar_bins_mass_temp_lum(1,solar_mass) == paired_stars(5, star_pair)
            T2 = solar_bins_mass_temp_lum(2,solar_mass);
            Bx2 = @(x) (1/(b_sigma*sqrt(2*pi)))*(exp(-0.5*((x-b_mu)/(b_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T2))-1));
            Vx2 = @(x) (1/(v_sigma*sqrt(2*pi)))*(exp(-0.5*((x-v_mu)/(v_sigma)).^2)).*(2.*hc.*c./x.^5).*(1./(e.^(hc./(x.*kb.*T2))-1));
            b_integral2 = integral(Bx2, 300, 1000);
            v_integral2 = integral(Vx2, 300, 1000);
            bv_ratio2 = b_integral2 / v_integral2;
            final_data_binary_stars(4,star_pair+1) = bvcoeff(1)+bvcoeff(2)*solar_bins_mass_temp_lum(2,solar_mass)^bvcoeff(3);
            final_data_binary_stars(7,star_pair+1) = solar_bins_mass_temp_lum(3,solar_mass);
        end
    end
    b_v1 = str2double(final_data_binary_stars(3,star_pair+1));
    b_v2 = str2double(final_data_binary_stars(4,star_pair+1));
    V1 = b_v1/(bv_ratio1-1);
    B1 = b_v1 + V1;
    V2 = b_v2/(bv_ratio2-1);
    B2 = b_v2 + V2;
%     final_data_binary_stars(5,star_pair+1) = (b_v1 + b_v2)/2;
    final_data_binary_stars(5,star_pair+1) = -2.5*log10((B1+B2)/(V1+V2));
    final_data_binary_stars(8,star_pair+1) = str2double(final_data_binary_stars(6,star_pair+1)) + str2double(final_data_binary_stars(7,star_pair+1));
end
figure
scatter(final_data_single_stars(2,:), final_data_single_stars(3,:),'x')
set(gca, 'YScale', 'log')
hold on
scatter(final_data_binary_stars(5,2:length(final_data_binary_stars)),final_data_binary_stars(8,2:length(final_data_binary_stars)),'.')

% b = [];
% for i = 2:length(final_data_binary_stars)
%     if final_data_binary_stars(1,i) > final_data_binary_stars(2,i)
%         b(i-1) = final_data_binary_stars(3,i);
%     end
%     if final_data_binary_stars(1,i) < final_data_binary_stars(2,i)
%         b(i-1) = final_data_binary_stars(4,i);
%     end
%     if final_data_binary_stars(1,i) == final_data_binary_stars(2,i)
%         b(i-1) = final_data_binary_stars(3,i);
%     end
% end
% scatter(b,final_data_binary_stars(8,2:length(final_data_binary_stars)),'d')

% Make finer bin sizes?