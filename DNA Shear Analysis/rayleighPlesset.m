function dR = rayleighPlesset(t,R)
format long
% t = 0.000001;
% R = [0.000001,0];
        %relayPlesset(1) = R'
        %relayPlesset(2) = R''
        %water @ 0-5 degrees Celsius
    rho = 996; %liquid density (kg/m3)
    sig = 0.075; %surface tension (N/m)
    mu = 0.798*10^-3; %liquid viscosity (Pa*s)
    P0 = 10^5; %atmospheric pressure (Pa)
    Pv = 4240; %vapor pressure (Pa)
    I = 20000; %intensity of ultrasound wave (W/m2)
    c = 1480; %speed of sound in water (m/s)
    freq = 22000; %frequency of ultrasound (Hz(1/s))
    A = sqrt((2*I)/((rho*c)*(freq^2)));
    Pamp = rho*c*A*freq;%pressure amplitude in ultrasound wave (Pa)
    %(2*I*rho*c)^(1/2)
    Pinf = P0*(1.0+A*sin(freq*t));%pressure of liquid surrounding bubble
    %P0-Pamp*cos(2*pi*freq*t); 
    R0 = 0.0001; %maximum radius of bubble (m)
    gamma = 1.3; %coefficient for adiabatic function for water steam
    Pi0 = Pinf; %pressure of saturated steam in bubble @ maximum radius
    Pi = Pv+Pi0*((R0/R(1))^(3*gamma)); %gas pressure inside bubble
    dR = zeros(2,1);
    dR(1) = R(2);
    dR(2) = -((3/(2*R(1)))*(R(2)^2))+((Pi-Pinf-((2*sig)/R(1)))/(rho*R(1)))-...
        ((4*mu*R(2))/(rho*(R(1)^2))); 
end