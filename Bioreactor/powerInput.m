Wvol = 0.02; %working volume, m^3
W = 0.018; %impellor width, m
Di = 0.215; %impellor diameter, m
Dt = 0.383; %vessel diameter, m
rho = 0.9821*1000; %density, kg/m3
mu = 0.0009; %viscosity, N s/m2
v = 0.9025*10^(-6); %kinematic viscosity, m2/s


N = 10.47; %rps
Re = ((Di^2)*(N))/v;

K4 = 1.1 + 4*(W/Dt) - 2.5*((Di/Dt - 0.5)^2) - 7*((W/Dt)^4);
K3 = 1.3 - 4*((W/Dt - 0.5)^2) - 1.14*(Di/Dt);
K2 = 10^K3;
K1 = 14 + (W/Dt)*(670*((Di/Dt - 0.6)^2) + 185);

NP = K1/Re + K2*(10^3) + 1.2*(Re^0.66)*(10^3) + 3.2*(Re^(0.66*K4));

NP_derived = 2.48*(Re^-0.31);

P = NP*(N^3)*(Di^5)*(rho); %Watts, kg*m2/s3

P_derived = NP_derived*(N^3)*(Di^5)*(rho);


%% 
%NAGATA TEST
rpm = [20; 40; 60; 80; 100; 120; 140; 160];
ReStore = [];
NPStore = [];
for i = 1:length(rpm)
    N = rpm(i)./60;
    Re = ((Di^2)*(N))/v;
    ReStore = [ReStore; Re];
    %Power Input
    K4 = 1.1 + 4*(W/Dt) - 2.5*((Di/Dt - 0.5)^2) - 7*((W/Dt)^4);
    K3 = 1.3 - 4*((W/Dt - 0.5)^2) - 1.14*(Di/Dt);
    K2 = 10^K3;
    K1 = 14 + (W/Dt)*(670*((Di/Dt - 0.6)^2) + 185);

    NP = K1/Re + K2*(10^3) + 1.2*(Re^0.66)*(10^3) + 3.2*(Re^(0.66*K4));
    NPStore = [NPStore; NP];
end