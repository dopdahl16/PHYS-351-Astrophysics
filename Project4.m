a_o = 269.25625604;
d_o = -18.9494472;
%Calibration_data
a = [269.3604167 269.2631917 269.2166375 269.1715333 269.0959833...
    269.0520083 269.4759750 269.2803917 269.2395667 269.2125750];
d = [-18.9406111 -18.9872000 -19.0009528 -19.1439694 -19.1496306...
    -19.0568000 -18.9058056 -19.0502278 -18.9541000 -19.0762750];
x_org = [746.58 547.57 462.38 230.97 112.45 153.63 957.99 500.73 550.72 369.82];
y_org = [459.64 493.95 523.01 347.2 420.52 614.49 389.36 376.04 571.94 409.03];
X = zeros(1,10);
Y = zeros(1,10);
for n = 1:10
    X(n) = -(cosd(d(n))*sind(a(n)-a_o))/(cosd(d_o)*cosd(d(n))*cosd(a(n) - a_o)...
        + sind(d_o)*sind(d(n)));
    Y(n) = -(sind(d_o)*cosd(d(n))*cosd(a(n)-a_o)-cosd(d_o)*sind(d(n)))/(cosd(d_o)*cosd(d(n))*cosd(a(n) - a_o)...
        + sind(d_o)*sind(d(n)));
end
plot(X, Y)
plot(x_org, y_org)
original = [x_org', y_org', ones(10,1)];
NEW = [X', Y'];
c_1 = zeros(1,8);
c_2 = zeros(1,8);
c_3 = zeros(1,8);
n_1 = zeros(1,8);
n_2 = zeros(1,8);
n_3 = zeros(1,8);
for n = 1:8
    A = inv(original(n:n+2,:))*NEW(n:n+2,:);
    %x transformation parameters
    c_1(n) = A(1,1);
    c_2(n) = A(2,1);
    c_3(n) = A(3,1);
    %y transformation parameters
    n_1(n) = A(1,2);
    n_2(n) = A(2,2);
    n_3(n) = A(3,2);
end
c_1 = mean(c_1);
c_2 = mean(c_2);
c_3 = mean(c_3);
n_1 = mean(n_1);
n_2 = mean(n_2);
n_3 = mean(n_3);
our_data = xlsread('Desktop\our stars');
X_trans = zeros(1, 2408);
a_dec = zeros(2408, 3);
for m = 1:2408
    x_orig = our_data(m,2);
    y_orig = our_data(m,3);
    X = x_orig*c_1 + y_orig*c_2 + c_3;
    Y = x_orig*n_1 + y_orig*n_2 + n_3;
    a_dec(m,1) = a_o + atand(-X/(cosd(d_o)-Y*sind(d_o)));
    a_dec(m,2) = asind((sind(d_o)+Y*cosd(d_o))/(sqrt(1 + X^2 + Y^2)));
    a_dec(m,3) = m;
end

% Retrieve data
[status,sheets] = xlsfinfo('M23refdata.xlsx');
M23RefImageData = xlsread('M23refdata.xlsx', char(sheets(1)));
WebDA = xlsread('M23refdata.xlsx', char(sheets(2)));

No_AUD = [M23RefImageData(:,1) M23RefImageData(:,4)];
No_VMag = [WebDA(:,2) WebDA(:,5)];

% Remove NaNs
row = 1;
while row <= length(No_AUD)
    if isnan(No_AUD(row,1)) || isnan(No_AUD(row,2))
        No_AUD(row,:) = [];
    else
        row = row + 1;
    end
end

row = 1;
while row <= length(No_VMag)
    if isnan(No_VMag(row,1)) || isnan(No_VMag(row,2))
        No_VMag(row,:) = [];
    else
        row = row + 1;
    end
end

% Cross Reference No_AUD and No_VMag
AUD_VMag = [];
for row = 1:length(No_AUD)
    ref_no = No_AUD(row,1);
    for line = 1:length(No_VMag)
        if No_VMag(line,1) == ref_no
            AUD_VMag = [AUD_VMag;ref_no No_AUD(row,2) No_VMag(line,2)];
        end
    end
end

% Graph VMag as a function of AUD
scatter(AUD_VMag(:,2),AUD_VMag(:,3))
hold on;

% Fit log to data
myfittype = fittype('a +b*log(x)',...
    'dependent', {'y'}, 'independent', {'x'},...
    'coefficients', {'a', 'b'});
myfit = fit(AUD_VMag(:,2),AUD_VMag(:,3),myfittype, 'StartPoint', [1,1]);
plot(myfit)

% Double check the fit is accurate
% hold on;
% coeffnames(myfit)
% coeffvalues(myfit)
% x = 1:100:1000000;
% f = 24.0495 + -1.1158*log(x);
% plot(x,f)