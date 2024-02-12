function [ Fps ] = forcepp( r,i,N,Nwall )
    global LX LX1 LX2 Ndtot Ntot
    global Nptot leff lp kBTp
%Initializing Array setup
    FS=zeros(2,Ntot);
    Fps=zeros(2,Ntot);
    Fpint=zeros(2,Ntot);
    diffrp=zeros(2,Ntot);
    absrp=zeros(1,Ntot);
    diffrvecp=zeros(2,Ntot);
    Fpint(1,i)=0;
    Fpint(2,i)=0;
    for j=i-1:2:i+1
        if i==Ndtot+1
            j=i+1;
        elseif i==Ndtot+Nptot
            j=i-1;
        end
        diffrp(1,i)=r(1,i)-r(1,j);
        if diffrp(1,i) > LX1
            diffrp(1,i)=LX - diffrp(1,i);
        elseif diffrp(1,i) < LX2
            diffrp(1,i)=abs(diffrp(1,i))-LX;
        end
        diffrp(2,i)=r(2,i)-r(2,j);
        absrp(i)=sqrt((diffrp(1,i))^2+(diffrp(2,i))^2);
        diffrvecp(1,i)=diffrp(1,i)./absrp(i);
        diffrvecp(2,i)=diffrp(2,i)./absrp(i);
%Spring Force between beads in a strand
        FS(1,i)=((-kBTp)/(4*leff))*(1-(absrp(i)/lp)^(-2)+(4*absrp(i)/lp)-...
        1)*diffrvecp(1,i);
        FS(2,i)=((-kBTp)/(4*leff))*(1-(absrp(i)/lp)^(-2)+(4*absrp(i)/lp)-...
        1)*diffrvecp(2,i);
        Fpint(1,i)=FS(1,i)+Fpint(1,i);
        Fpint(2,i)=FS(2,i)+Fpint(2,i);
        if j==i+1
            break
        elseif i==Ndtot+Nptot
            break
        end
    end
    Fps(1,i)=Fpint(1,i);
    Fps(2,i)=Fpint(2,i);
end