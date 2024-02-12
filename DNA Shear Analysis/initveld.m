function [ v ] = initveld( N,Nwall,Np,Npb )
    global kBT md Nptot Ndtot Ntot
    global velmxd
%Initiate Arrays
    v=zeros(2,Ntot);
    for i=1:1:N
        v(1,i)=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
        v(2,i)=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
% Velocity check to ensure velocity does not exceed max velocity velmxd
        vavg=(v(1,i))^2+(v(2,i))^2;
        if vavg > velmxd
            vavg=sqrt(velmxd/vavg);
            v(1,i)=v(1,i)*vavg;
            v(2,i)=v(2,i)*vavg;
        end
    end
% To set Total Momentum equals Zero
    momxd=0;
    momyd=0;
    for i=1:1:N
        momxd=momxd+v(1,i);
        momyd=momyd+v(2,i);
    end
        momxd=momxd/N;
        momyd=momyd/N;
    for i=1:1:N
        v(1,i)=v(1,i)-momxd;
        v(2,i)=v(2,i)-momyd;
    end
%%%%%%%%%%%%%%%%%%%
% Wall Velocity
%%%%%%%%%%%%%%%%%%%
    for i=N+1:1:Ndtot
        v(1,i)=0;
        v(2,i)=0;
    end
%%%%%%%%%%%%%%%%%%%
% DNA Velocity
%%%%%%%%%%%%%%%%%%%
    for i=Ndtot+1:1:Ndtot+Nptot
        v(1,i)=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
        v(2,i)=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
% Velocity check to ensure velocity does not exceed max velocity velmxd
        vavg=(v(1,i))^2+(v(2,i))^2;
        if vavg > velmxd
            vavg=sqrt(velmxd/vavg);
            v(1,i)=v(1,i)*vavg;
            v(2,i)=v(2,i)*vavg;
        end
    end