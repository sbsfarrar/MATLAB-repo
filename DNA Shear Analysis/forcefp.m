function [Fp ] = forcefp(r,v,i,N,Nwall)
    global LX LX1 LX2
    global kBT rc rc2 s
    global Ndtot Nptot Ntot
    global sigma gamma g delt
    global afp apw
%Initiate Arrays
    diffr=zeros(2,Ntot); diffv=zeros(2,Ntot);
    absr=zeros(1,Ntot); absr2=zeros(1,Ntot); absv=zeros(1,Ntot);
    diffrvec=zeros(2,Ntot); diffvvec=zeros(2,Ntot);
    FCon=zeros(2,Ntot); FDis=zeros(2,Ntot);
    FRan=zeros(2,Ntot);
    Fp=zeros(2,Ntot);
    dotrv=zeros(2,Ntot);
    Fint(1,i)=0;
    Fint(2,i)=0;
    Fintw(1,i)=0;
    Fintw(2,i)=0;
    Fintp(1,i)=0;
    Fintp(2,i)=0;
%%%%%%%%%%%%%%%%%%%%
% Fluid Particles
%%%%%%%%%%%%%%%%%%%%
    for j=1:1:N
%Distance between two particles with x and v components
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1
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
        FCon(1,i)=afp*(1-absr(i)).*diffrvec(1,i);
        FCon(2,i)=afp*(1-absr(i)).*diffrvec(2,i);
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
        dotrv(i)=(diffrvec(1,i).*diffv(1,i))+(diffrvec(2,i).*diffv(2,i));
        FDis(1,i)=-gamma*wD*dotrv(i).*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dotrv(i).*diffrvec(2,i);
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
        Fint(1,i)=Fint(1,i)+(FCon(1,i)+FDis(1,i)+FRan(1,i));
        Fint(2,i)=Fint(2,i)+(FCon(2,i)+FDis(2,i)+FRan(2,i));
    end
%%%%%%%%%%%%%%%%%%%%%%%
% Wall Particles
%%%%%%%%%%%%%%%%%%%%%%%
    for j=N+1:1:N+Nwall
%Distance between two particles with x and v components
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1
            diffr(1,i)= diffr(1,i)-LX;
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
        FCon(1,i)=apw*(1-absr(i)).*diffrvec(1,i);
        FCon(2,i)=apw*(1-absr(i)).*diffrvec(2,i);
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
        dotrv(i)=(diffrvec(1,i).*diffv(1,i))+(diffrvec(2,i).*diffv(2,i));
        FDis(1,i)=-gamma*wD*dotrv(i).*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dotrv(i).*diffrvec(2,i);
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
        Fintw(1,i)=Fintw(1,i)+(FCon(1,i)+FDis(1,i)+FRan(1,i));
        Fintw(2,i)=Fintw(2,i)+(FCon(2,i)+FDis(2,i)+FRan(2,i));
    end
%%%%%%%%%%%%%%%%%%%%%%
% DNA Particles
%%%%%%%%%%%%%%%%%%%%%%
    for j=Ndtot+1:1:Ndtot+Nptot
        if j==i
            continue
        end
%Distance between two particles with x and v components
        diffr(1,i)=r(1,i)-r(1,j);
        if diffr(1,i) > LX1
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
        absr(i)=sqrt((diffr(1,i))^2+(diffr(2,i))^2);
        diffrvec(1,i)=diffr(1,i)./absr(i);
        diffrvec(2,i)=diffr(2,i)./absr(i);
%Velocity between two particles with x and v components
        diffv(1,i)=v(1,i)-v(1,j);
        diffv(2,i)=v(2,i)-v(2,j);
        absv(i)=sqrt((diffv(1,i))^2+(diffv(2,i))^2);
        diffvvec(1,i)=diffv(1,i)./absv(i);
        diffvvec(2,i)=diffv(2,i)./absv(i);
%Conservative Force- Repulsive Force
        if abs(i-j)>4
            app=2;
        else
            app=0;
        end
        FCon(1,i)=app*(1-absr(i))*diffrvec(1,i);
        FCon(2,i)=app*(1-absr(i))*diffrvec(2,i);
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
        FDis(1,i)=-gamma*wD*dot(diffrvec(1,i),diffv(1,i))*diffrvec(1,i);
        FDis(2,i)=-gamma*wD*dot(diffrvec(2,i),diffv(2,i))*diffrvec(2,i);
        FRan(1,i)=sigma*wR*theta*diffrvec(1,i);
        FRan(2,i)=sigma*wR*theta*diffrvec(2,i);
        Fintp(1,i)=Fintp(1,i)+FCon(1,i)+FDis(1,i)+FRan(1,i)*delt^(-0.5);
        Fintp(2,i)=Fintp(2,i)+FCon(2,i)+FDis(2,i)+FRan(2,i)*delt^(-0.5);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Forces on Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fext(1,i)=g;
    Fext(2,i)=0;
    Fp(1,i)=Fint(1,i)+Fext(1,i)+Fintw(1,i)+Fintp(1,i);
    Fp(2,i)=Fint(2,i)+Fext(2,i)+Fintw(2,i)+Fintp(2,i);
end