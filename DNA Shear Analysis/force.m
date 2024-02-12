function [ F ] = force( r,v,i,N,Nwall )
    global kBT aff afw afp rc rc2 s delt
    global LX LX1 LX2
    global Ndtot Nptot Ntot
    global sigma gamma g
%Initiate Arrays
    diffr=zeros(2,N); diffv=zeros(2,N);
    absr=zeros(1,N); absr2=zeros(1,N); absv=zeros(1,N);
    diffrvec=zeros(2,N); diffvvec=zeros(2,N);
    FCon=zeros(2,N); FDis=zeros(2,N);
    FRan=zeros(2,N);
    F=zeros(2,Ntot);
    dotrv=zeros(2,N);
    Fint=zeros(2,N); Fintw=zeros(2,N); Fintp=zeros(2,N);
    Fext(1,1:N)=g;
    Fext(2,1:N)=0;
    Fint(1,i)=0;
    Fint(2,i)=0;
    Fintw(1,i)=0;
    Fintw(2,i)=0;
    Fintp(1,i)=0;
    Fintp(2,i)=0;
    for j=1:1:N
        if j==i
            continue
        end
%Distance between two particles with x and v components
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1 % Periodic Boundary Conditions
            diffr(1,i)=diffr(1,i)-LX;
        elseif diffr(1,i) < LX2
            diffr(1,i)=LX - abs(diffr(1,i));
        end
        diffr(2,i)=r(2,i)-r(2,j);
        if abs(diffr(1,i))>rc
            continue
        end
        if abs(diffr(2,i))>rc
            continue
        end
        absr(i)=sqrt((diffr(1,i)).^2+(diffr(2,i)).^2);
        absr2(i)=(absr(i)).^2;
        if absr2(i)>rc2
            continue
        end
        diffrvec(1,i)=diffr(1,i)./absr(i);
        diffrvec(2,i)=diffr(2,i)./absr(i);
%Velocity between two particles with x and v components
        diffv(1,i)=v(1,i)-v(1,j);
        diffv(2,i)=v(2,i)-v(2,j);
        absv(i)=sqrt((diffv(1,i)).^2+(diffv(2,i)).^2);
        diffvvec(1,i)=diffv(1,i)./absv(i);
        diffvvec(2,i)=diffv(2,i)./absv(i);
%Conservative Force- Repulsive Force
        FCon(1,i)=aff*(1-absr(i)).*diffrvec(1,i);
        FCon(2,i)=aff*(1-absr(i)).*diffrvec(2,i);
        if abs(absr(i))<=rc
            wD=(1-absr(i)/rc)^s;
        else
            wD=0;
        end
        wR=sqrt(wD);
        theta= sqrt((-2)*log(rand))*cos(2*pi*rand);
        if theta > 6
            theta=sign(theta)*6;
        end
%Dissipative Force
        dotrv(i)=(diffrvec(1,i).*diffv(1,i))+(diffrvec(2,i).*diffv(2,i));
        FDis(1,i)=-gamma*wD*dotrv(i).*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dotrv(i).*diffrvec(2,i);
%Random Force
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
%Internal Forces
        Fint(1,i)=Fint(1,i)+(FCon(1,i)+FDis(1,i)+FRan(1,i));
        Fint(2,i)=Fint(2,i)+(FCon(2,i)+FDis(2,i)+FRan(2,i));
    end
%%%%%%%%%%%%%%%%%%%%%%%
% Wall Particles
%%%%%%%%%%%%%%%%%%%%%%%
    for j=N+1:1:N+Nwall
%Distance between two particles with x and v components
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1 % Periodic Boundary Conditions
            diffr(1,i)=diffr(1,i)-LX;
        elseif diffr(1,i) < LX2
            diffr(1,i)=LX - abs(diffr(1,i));
        end
        diffr(2,i)=r(2,i)-r(2,j);
        if abs(diffr(1,i))>rc
            continue
        end
        if abs(diffr(2,i))>rc
            continue
        end
        absr(i)=sqrt((diffr(1,i)).^2+(diffr(2,i)).^2);
        absr2(i)=(absr(i)).^2;
        if absr2(i)>rc2
            continue
        end
        diffrvec(1,i)=diffr(1,i)./absr(i);
        diffrvec(2,i)=diffr(2,i)./absr(i);
%Velocity between two particles with x and v components
        diffv(1,i)=v(1,i)-v(1,j);
        diffv(2,i)=v(2,i)-v(2,j);
        absv(i)=sqrt((diffv(1,i)).^2+(diffv(2,i)).^2);
        diffvvec(1,i)=diffv(1,i)./absv(i);
        diffvvec(2,i)=diffv(2,i)./absv(i);
%Conservative Force- Repulsive Force
        FCon(1,i)=afw*(1-absr(i)).*diffrvec(1,i);
        FCon(2,i)=afw*(1-absr(i)).*diffrvec(2,i);
        if abs(absr(i))<=rc
            wD=(1-absr(i)/rc)^s;
        else
            wD=0;
        end
%Weight Functions and Coefficients of FD and FR
        wR=sqrt(wD);
        theta= sqrt((-2)*log(rand))*cos(2*pi*rand);
        if theta > 6
            theta=sign(theta)*6;
        end
%Dissipative Force
        dotrv(i)=(diffrvec(1,i).*diffv(1,i))+(diffrvec(2,i).*diffv(2,i));
        FDis(1,i)=-gamma*wD*dotrv(i).*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dotrv(i).*diffrvec(2,i);
%Random Force
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
%Internal Forces
        Fintw(1,i)=Fintw(1,i)+(FCon(1,i)+FDis(1,i)+FRan(1,i));
        Fintw(2,i)=Fintw(2,i)+(FCon(2,i)+FDis(2,i)+FRan(2,i));
    end
%%%%%%%%%%%%%%%%%%%%%%
% DNA Particles
%%%%%%%%%%%%%%%%%%%%%%
    for j=Ndtot+1:1:Ndtot+Nptot
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1 % Periodic Boundary Conditions
            diffr(1,i)=diffr(1,i)-LX;
        elseif diffr(1,i) < LX2
            diffr(1,i)=LX - abs(diffr(1,i));
        end
        diffr(2,i)=r(2,i)-r(2,j);
        if abs(diffr(1,i))>rc % Setting neighboring particles
            continue
        end
        if abs(diffr(2,i))>rc
            continue
        end
        absr(i)=sqrt((diffr(1,i)).^2+(diffr(2,i)).^2);
        absr2(i)=(absr(i)).^2;
        if absr2(i)>rc2
            continue
        end
        diffrvec(1,i)=diffr(1,i)./absr(i);
        diffrvec(2,i)=diffr(2,i)./absr(i);
%Velocity between two particles with x and v components
        diffv(1,i)=v(1,i)-v(1,j);
        diffv(2,i)=v(2,i)-v(2,j);
        absv(i)=sqrt((diffv(1,i)).^2+(diffv(2,i)).^2);
        diffvvec(1,i)=diffv(1,i)./absv(i);
        diffvvec(2,i)=diffv(2,i)./absv(i);
%Conservative Force- Repulsive Force
        FCon(1,i)=afp*(1-absr(i)).*diffrvec(1,i);
        FCon(2,i)=afp*(1-absr(i)).*diffrvec(2,i);
        if abs(absr(i))<=rc
            wD=(1-absr(i)/rc)^s;
        else
            wD=0;
        end
        wR=sqrt(wD);
        gamma=(sigma^2)/(2*kBT);
        theta= sqrt((-2)*log(rand))*cos(2*pi*rand);
        if theta > 6
            theta=sign(theta)*6;
        end
%Dissipative Force
        dotrv(i)=(diffrvec(1,i).*diffv(1,i))+(diffrvec(2,i).*diffv(2,i));
        FDis(1,i)=-gamma*wD*dotrv(i).*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dotrv(i).*diffrvec(2,i);
%Random Force
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
%Internal Polymer Forces
        Fintp(1,i)=Fintp(1,i)+FCon(1,i)+FDis(1,i)+FRan(1,i)*delt^(-0.5);
        Fintp(2,i)=Fintp(2,i)+FCon(2,i)+FDis(2,i)+FRan(2,i)*delt^(-0.5);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Forces on Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F(1,i)=Fint(1,i)+Fext(1,i)+Fintw(1,i)+Fintp(1,i);
    F(2,i)=Fint(2,i)+Fext(2,i)+Fintw(2,i)+Fintp(2,i);
end