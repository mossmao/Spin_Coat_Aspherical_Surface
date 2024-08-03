%% define coating parameters
rpm = 1000;    
w = rpm * 2* pi / 60; % rad/s
g = 9.8;  % m/s^2  

rho = 1090; % kg/m^3, resin density
eta = 0.370;  % Pa.s, resin viscosity

tend = 20; % seconds, total coating time
h0=0.25e-3; % meter, the initial thickness at the center

%% define an asperhic lens profile
R = 15*1e-3; % meter
conic = 0;          % conic
A = [0, 0, 0, 0];   % coefficents

% descretize the profile into points
N = 1e6;
r = linspace(0,R*sqrt(1/(1+conic)*0.6),N);
z = -r.^2./(R*(1+sqrt(1-(1+conic)*r.^2/R^2)));

for ii = 1:length(A)
    z = z - A(ii)*r.^(2*ii+2);
end

%% compute the profile
%predict the coating thickness of the initial design
[h_m, ~,~, theta] = predict_coat_thickness(r,z,eta,rho,w,g,h0,tend);

% revise the design to compensate for the coating thickness
r_corrected = r - h_m.*sin(theta);
z_corrected = z - h_m.*cos(theta);

z_new = interp1(r_corrected,z_corrected,r,'linear','extrap');

% predict the coated profile of corrected substrated 
[h_m, r_coated,z_coated] = predict_coat_thickness(r_corrected,z_corrected,eta,rho,w,g,h0,tend);
z_new_coated = interp1(r_coated,z_coated,r,'linear','extrap');


figure; hold on;
set(gcf,'Units','Inches')
set(gcf, 'Position', 2+[0 0 2.5 2]);
set(gcf,'Resize','off');
set(gca,'FontSize',14)
set(gca,'FontName','Times New Roman')

plot(r*1000,z*1000,r*1000,z_new*1000,r*1000,z_new_coated*1000);
legend('Designed','Compensated','Coated');
xlabel('r (mm)');
ylabel('z (mm)');
r_mm = r*1000;
h_mm = h_m*1000;

raw_data = [r*1000;z*1000;r*1000;z_new*1000;r*1000;z_new_coated*1000];
raw_data = raw_data(:,1:100:end);
raw_data

figure; hold on;
set(gcf,'Units','Inches')
set(gcf, 'Position', 2+[0 0 2.5 2]);
set(gcf,'Resize','off');
set(gca,'FontSize',14)
set(gca,'FontName','Times New Roman')

yyaxis right;
plot(r*1000,(z_new_coated-z)*10^6);

our = [r*1000;(z_new_coated-z)*10^6];
ylim([min((z_new_coated-z)*10^6),0.1]);
ylabel('Our dimension error (\mum)');
xlabel('r (mm)');
mean(abs(z_new_coated-z))
%% uniform compensation
% revise the design to compensate for the coating thickness
r_corrected = r - h_m(1).*sin(theta);
z_corrected = z - h_m(1).*cos(theta);

z_new = interp1(r_corrected,z_corrected,r,'linear','extrap');

% predict the coated profile of corrected substrated 
[h_m, r_coated,z_coated] = predict_coat_thickness(r_corrected,z_corrected,eta,rho,w,g,h0,tend);
z_new_coated = interp1(r_coated,z_coated,r,'linear','extrap');


yyaxis left;
plot(r*1000,(z_new_coated-z)*10^6);
mean(abs(z_new_coated-z))
legend('Uniform','Our','location','nw');
ylabel('Conventional dimension error (\mum)');
old = [r*1000;(z_new_coated-z)*10^6];
raw_data2 = [our;old];
raw_data2 = raw_data2(:,1:100:end);
raw_data2
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [2.5 2]);
% set(gca,'FontSize',10)
% set(gcf,'Units','Inches')

function [h_m, r_coated,z_coated, theta] = predict_coat_thickness(r_m,z_m,eta,rho,w,g,h0,tend)
N = length(r_m);

% estimate local Radius at the center
r_L = r_m(1:1000)';
z_L = z_m(1:1000)';

R = (2*(z_L(1)-z_L)) \ (r_L.^2+(z_L(1)-z_L).^2);

tc = 0.75*eta/h0^2/(rho*(w*w+g/R));

r = r_m;
z = z_m;
dr = diff(r);
dz = diff(z);

ds = sqrt(dr.^2+dz.^2);

s = cumsum(ds);
theta = atan2(-dz,dr);
r = r(2:end);
z = z(2:end);

fs = rho/eta*(w^2*r.^2.*cos(theta)+r.*g.*sin(theta));
us = fs.^(-1/3).*(cumsum(r.*fs.^(-1/3).*ds)).^(1/2);

us(1:(N/100)) = us(N/100); % this line is to remove the singlarity at initial points
fun_t = 1./sqrt(tend+tc);

h_m = fun_t*us;
h_m = [h_m(1) h_m];

% offset the z
theta = [theta(1) theta];

r_coated = r_m + h_m.*sin(theta);
z_coated = z_m + h_m.*cos(theta);
end