%% Relay-Plesset Equation Solved via Runge-Kutta
    %RUNGE-KUTTA-FEHLBERG
    % It calculates ODE using Runge-Kutta 4th and 5th order method
    % from H.Alehossein and Z. Qin 2007
h=2*10^-9; % step size (s)
t = 0.00000:h:0.00008;  % span(s)
R = zeros(2,length(t)); %R(1) = R, R(2) = R'
E = zeros(2,length(t)); %Error
R(1,1) = 0.00001;    %initial R @ t0 (m)                                
R(2,1) = 0;     %initial R'@ t0
                  
for counti=1:(length(t)-1)                                   % calculation loop
    k_1 = h*rayleighPlesset(t(counti),R(:,counti));
    k_2 = h*rayleighPlesset((t(counti)+(h/4)),(R(:,counti)+(k_1/4)));
    k_3 = h*rayleighPlesset((t(counti)+((3*h)/8)),(R(:,counti)+((3*k_1)/32)+((9*k_2)/32)));
    k_4 = h*rayleighPlesset((t(counti)+((12*h)/13)),(R(:,counti)+((1932*k_1)/2197)-((7200*k_2)/2197)+...
        ((7296*k_3)/2197)));
    k_5 = h*rayleighPlesset((t(counti)+h),(R(:,counti)+((439*k_1)/216)-8*k_2+((3680*k_3)/513)-...
        ((845*k_4)/4104))); 
    k_6 = h*rayleighPlesset((t(counti)+(h/2)),(R(:,counti)-((8*k_1)/27)+2*k_2-((3544*k_3)/2565)+...
        ((1859*k_4)/4104)-((11*k_5)/40))); 

    %R(:,i+1) =
    %R(:,i)+((25*k_1)/216)+((1408*k_3)/2565)+((2197*k_4)/4101)-(k_5/5);
    %main equation RKF of order 4
        
    R(:,counti+1) = R(:,counti)+(((16*k_1)/135)+((6656*k_3)/12825)+((28561*k_4)/56430)-...
        ((9*k_5)/50)+((2*k_6)/55));  % main equation RKF of order 5
    E(:,counti) = (k_1/360)-((128*k_3)/4275)-((2197*k_4)/75240)+(k_5/50)+((2*k_6)/55);%error approx
end
%adiabatic approximation: 
% syms Pi0 R0 gamma
% Pi = Pi0*((R0/R)^(3*gamma));
%% Absolute Value of Radial Speed in Liquid Jet Stream
    %V = radial speed of liquid flux at distance r from the center of the 
    %    bubble (r>R);
r = linspace(0.0000001,0.000003,length(R));
G = zeros(length(r),length(r));  %|dV/dr|
T = zeros(length(r),length(r));
rr = zeros(length(r),length(r));
% [rr,tt] = meshgrid(r,t);
% GG = 2*((R1.^2)./(rr.^3)).*abs(R2);
for counti = 1:length(r)
    for countj = 1:length(r)
        G(counti,countj) = 2*((R(1,countj)^2)/(r(counti)^3))*abs(R(2,countj));
        T(counti,countj) = t(countj);
        rr(counti,countj) = r(counti);
    end
end
G = G(1:50:numel(r),1:50:numel(r));
T = T(1:50:numel(r),1:50:numel(r));
rr = rr(1:50:numel(r),1:50:numel(r));
mesh(rr,T,G);
xlabel('rr');
ylabel('T');
zlabel('G');
% %%
% %Stockes' Law (3):
% syms F Rd
% F = 6*pi*mu*Rd*V;
% %ith chain model (4):
% syms m k l0 xiplus xi ximinus alpha Va x
% alpha =6*pi*mu*Rd;
% ode2 = m*diff(xi,t,2) == ((k*(xiplus - xi - l0))/l0)-((k*(xi-ximinus-l0))/l0)+...
%     alpha*(Va - diff(x));
% 
% 
% %%
% %speed of mechanochemical reaction (5):
% syms kp T h delS delEf
% kp = ((k*T)/(2*pi*h))*exp(delS/k)*exp(-(delEf)/(k*T)); 
% %energetic barrier for formation of the activated complex (6):
% syms Estar E0 f delQ
% delEf = Estar - E0 - f*delQ;
% %speed of reaction in presence of stretching force, F (7):
% syms kf k0
% kf = k0*exp((f*delQ)/(k*T));
% %Morse potential describing deformation of a covalent bond (8): 
% syms Vm delEstar
% Vm = delEstar; %Vm(delQ) = D*((1-exp(-beta*delQ))^2)
% %stretching of a specific internucleotide bond (9): 
% syms Fstretch
% f = gamma*Fstretch;
% %portion of cleaved DNA fragments at a particular stage of bubble
% %implosion (10): 
% syms krel n0 beta 
% krel = int(int(4*pi*(r^2)*n0*beta*k0*exp((gamma*Fstretch)/k*T),r),t)/...
%     int(4*pi*(r^2)*n0,r);
% %kN and kS approximations (11): 