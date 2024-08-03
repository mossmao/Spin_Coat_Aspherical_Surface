%% define coating parameters
rpm = 1000;    
w = rpm * 2* pi / 60; % rad/s
g = 9.8;  % m/s^2  

rho = 1090; % kg/m^3, resin density
eta = 0.370;  % Pa.s, resin viscosity

tend = 100; % seconds, total coating time
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

%% plot the results
tc = 0.75*eta/h0^2/(rho*(w*w+g/R));
% plot the thickness at the center that changes with time
figure;
t = linspace(0,tend,1000);
thickness_at_center = h0*sqrt(tc)./sqrt(t+tc)*1000;

plot(t, thickness_at_center);
xlabel('Time (s)');
ylabel('Coating thickness at center(mm)');



% plot the thickness profile at the end coating time
figure;
t = tend;
fun_t = 1./sqrt(t+tc)*1000;

thickness_at_tend = fun_t*us;
plot(r*1000, thickness_at_tend);
xlabel('Radius r (mm)');
ylabel('Coating thickness at t_{end}(mm)');