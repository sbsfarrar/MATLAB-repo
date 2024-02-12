tic
global kBT aff aww afw rc rc2 rcw s delt
global LX LY LX1 LX2 LY1 LY2
global Nwall Nwall2 Q Qwall Ndtot Nptot Ntot
global velmxd
global A Qm1 nx ny nwallx nwally
global sigma gamma rho md g lambda
global afp apw Np Npb lseg nyp leff L lp kBTp
% PARAMETERS
%DPD constants
dm=2; %Dimensions
kBT=1; %=kB*Temp
sigma=3;
gamma=4.5;
lambda=0.65;
rho=4;
aff=(75*kBT)/rho; % =18.75
aww=5.0; % =5
afw=sqrt(aff*aww); % =9.682
rc=1;
rc2=rc^2;
s=2;
rl=1.5; % Verlet Neighbour List Method r<rc<=rl
delt=0.02; %time step
tf=30; %Number of time steps (t>1350)
ti=0.; %inital time
g=10; %Driving force in x direction
M=1; % Mass density of DPD particles
md=1;
velmxd=sqrt(2*kBT/md); % Maximum velocity of particles
%-------------------------------------------------------------
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPD Fluid Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1458; % Number of particles
LX=12;
LY=30;
rcw=0.005*LY;
ndensdn=1; % Non-dimensional Number density of DPD
dc=0.4; % Diameter of Particle
ndensd=N/(LX*LY); % Number density of dpd particles
vdensd=ndensdn*pi/4; % Volumetric Fraction
Q=2*sqrt(N/2);
A=sqrt(1/ndensd); % Number density
Qm1=(Q-1);
nx=(LX)/Q;
ny=(LY)/Q;
LX1=LX/2;
LX2=-(LX/2);
LY1=LY/2;
LY2=(-LY/2);
bc1=(LY1)-rcw;
bc2=(LY2)+rcw;
bins=50;
tbins=bins+2; %Total Bins
nybins=(LY-(2*rcw))/bins;
nymat=[(LY2+rcw/2) (LY2+rcw+(nybins/2):nybins:LY1-rcw) (LY1-rcw/2)];
bnbtm=rcw-(LY1); %Bottom bin near wall ny=0.005*LY
bntp=(LY1)-rcw; %Top bin near wall ny=0.005*LY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wall Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nwall=400; %Number of Wall particles
Nwall2=Nwall/2;
Qwall=sqrt(Nwall);
nwallx=LX/Nwall2;
nwally=LY/Nwall2;
Ndtot=Nwall+N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNA Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np=30; % Number of strands
Npb=10; %Number of beads in the strand=81
leff=0.053; %effective length;
kBTp=1; %uJ kBTp=4.115*10^(-14)erg Erg=1*10^-7 J
L=67.2; %total length of the DNA strand 67.2 um
lp=0.81; % for Npb=81 and lp=L/(Npb-1);
lseg=0.4; % Initial distance between beads
mp=0.25; % g/cm3
viscp=2.588;
nvp=1.235; %cP um/s
Nptot=Np*Npb; % Total number of beads in all strands
Ntot=Ndtot+Nptot;
afp=2; % repulsion force between fluid and polymer
apw=sqrt(2*aww);
nyp=LY/Np;
%Initiate Arrays
newv=zeros(2,Ntot);
newF=zeros(2,Ntot); predv=zeros(2,Ntot);
vavgd2=zeros(1,tbins); vavgd=zeros(1,tbins);
vnum=zeros(1,tbins);
diffrpl=zeros(2,Nptot);
absrpl=zeros(1,Nptot);
constl=zeros(1,Nptot);
vecl=zeros(2,Nptot);
%Initiate Conditions
r=initposd(N,Nwall);
newr=initposd(N,Nwall);
v=initveld(N,Nwall,Np,Npb);
vnum(1,1:tbins)=0;
vavgd(1,1:tbins)=0;
F=0; Fp=0; Fps=0;
%Initial Force
for i=1:1:N
    F=F+force(r,v,i,N,Nwall);
end
for i=Ndtot+1:1:Ndtot+Nptot
    Fp=Fp+forcefp(r,v,i,N,Nwall);
    Fps=Fps+forcepp(r,i,N,Nwall);
end
F=F+Fp+Fps;
%Modified Velocity Verlet
for t=ti:delt:tf
    t;
    for i=1:1:N
        for k=1:2
            if abs(r(2,i))<bc1
                newr(k,i)=r(k,i)+delt*v(k,i)+(1/2)*delt^2*(1/M)*F(k,i);
                predv(k,i)=v(k,i)+ lambda*delt*(1/M)*F(k,i);
            else
                newr(k,i)=r(k,i)+delt*v(k,i);
                predv(k,i)=v(k,i);
            end
        end
%Setting Periodic Boundary Conditions
        if newr(1,i)>= LX1
            newr(1,i)=newr(1,i)-LX;
        elseif newr(1,i)<= LX2
            newr(1,i)=newr(1,i)+LX;
        end
    end
% When particles are close to the wall particles
    for i=1:1:N
        if newr(2,i)>=bc1 % Top wall
            n=-1;
            vRx=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            vRy=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            if abs(v(1,i))>=abs(vRx)
                newv(1,i)=vRx;
            else
                newv(1,i)=v(1,i);
            end
            if abs(v(2,i))>=abs(vRy)
                newv(2,i)=vRy+n*(sqrt((n*vRy)^2)-(n*vRy));
                if newv(2,i)>0
                    disp('positive at upper wall')
                    newv(2,i)=newv(2,i)*(-1);
                end
            else
                newv(2,i)=n*abs(v(2,i));
            end
            newF=zeros(2,Ntot);
            newr(2,i)=bc1;
        elseif newr(2,i)<=bc2 % Bottom wall
            n=1;
            vRx=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            vRy=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            if abs(v(1,i))>=abs(vRx)
                newv(1,i)=vRx;
            else
                newv(1,i)=v(1,i);
            end
            if abs(v(2,i))>=abs(vRy)
                newv(2,i)=vRy+n*(sqrt((n*vRy)^2)-(n*vRy));
                if newv(2,i)<0
                    disp('negative at lower wall')
                    newv(2,i)=newv(2,i)*(-1);
                end
            else
                newv(2,i)=n*abs(v(2,i));
            end
            newF=zeros(2,Ntot);
            newr(2,i)=bc2;
        else
            newF=force(newr,predv,i,N,Nwall);
            newv(1,i)=v(1,i)+(1/2)*delt*(1/M)*(F(1,i)+newF(1,i));
            newv(2,i)=v(2,i)+(1/2)*delt*(1/M)*(F(2,i)+newF(2,i));
        end
% Velocity check to ensure velocity does not exceed
% maximum velocity velmxd
        vavg=(newv(1,i))^2+(newv(2,i))^2;
        if vavg > velmxd^2
            vavg2=sqrt(velmxd^2/vavg);
            newv(1,i)=newv(1,i)*vavg2;
            newv(2,i)=newv(2,i)*vavg2;
        end
        r(1,i)=newr(1,i);
        r(2,i)=newr(2,i);
        v(1,i)=newv(1,i);
        v(2,i)=newv(2,i);
        F(1,i)=newF(1,i);
        F(2,i)=newF(2,i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNA Particle interaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=Ndtot+1:1:Ndtot+Nptot
        for k=1:2
            if abs(r(2,i))<bc1
                newr(k,i)=r(k,i)+delt*v(k,i)+(1/2)*delt^2*(1/M)*F(k,i);
                predv(k,i)=v(k,i)+ lambda*delt*(1/M)*F(k,i);
            else
                newr(k,i)=r(k,i)+delt*v(k,i);
                predv(k,i)=v(k,i);
            end
        end
% Left beads are recorded and the right beads are adjusted with
% length less than or equal to 0.01
        diffrpl(1,i)=newr(1,i)-newr(1,i-1);
        if diffrpl(1,i) > LX1
            diffrpl(1,i)=diffrpl(1,i)-LX;
        elseif diffrpl(1,i) < LX2
            diffrpl(1,i)=LX - abs(diffrpl(1,i));
        end
        diffrpl(2,i)=newr(2,i)-newr(2,i-1);
        absrpl(i)=sqrt((diffrpl(1,i)).^2+(diffrpl(2,i)).^2);
        constl(i)=absrpl(i)/lp;
        vecl(1,i)=diffrpl(1,i)./absrpl(i);
        vecl(2,i)=diffrpl(2,i)./absrpl(i);
% Constraint to prevent length of segments near fixed end exceeding lp
        if i~=Ndtot+1
            if abs(constl(i)) > 0.95
                newr(1,i)=newr(1,i-1)+lp*vecl(1,i);
                newr(2,i)=newr(2,i-1)+lp*vecl(2,i);
            end
        end
%Setting Periodic Boundary Conditions
        if newr(1,i)>= LX1
            newr(1,i)=newr(1,i)-LX;
        elseif newr(1,i)<= LX2
            newr(1,i)=newr(1,i)+LX;
        end
    end
    for i=Ndtot+1:1:Ndtot+Nptot
        if newr(2,i)>=bc1 % Top wall
            n=-1;
            vRx=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            vRy=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            newv(1,i)=vRx;
            newv(2,i)=vRy+n*(sqrt((n*vRy)^2)-(n*vRy));
            if newv(2,i)>0
                disp('positive at upper wall')
                newv(2,i)=newv(2,i)*(-1);
            end
            newF=zeros(2,Ntot); %function Force
            newr(2,i)=0;
        elseif newr(2,i)<=bc2 % Bottom wall
            n=1;
            vRx=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            vRy=sqrt((-2)*(kBT/md)*log(rand))*cos(2*pi*rand);
            newv(1,i)=vRx;
            newv(2,i)=vRy+n*(sqrt((n*vRy)^2)-(n*vRy));
            if newv(2,i)<0
                disp('negative at lower wall')
                newv(2,i)=newv(2,i)*(-1);
            end
            newF=zeros(2,Ntot); %function Force
            newr(2,i)=0;
        else
            newF=forcefp(newr,predv,i,N,Nwall)+forcepp(newr,i,N,Nwall);
            newv(1,i)=v(1,i)+(1/2)*delt*(1/M)*(F(1,i)+newF(1,i));
            newv(2,i)=v(2,i)+(1/2)*delt*(1/M)*(F(2,i)+newF(2,i));
        end
        vavg=(newv(1,i))^2+(newv(2,i))^2;
        if vavg > velmxd^2
            vavg2=sqrt(velmxd^2/vavg);
            newv(1,i)=newv(1,i)*vavg2;
            newv(2,i)=newv(2,i)*vavg2;
        end
        r(1,i)=newr(1,i);
        r(2,i)=newr(2,i);
        v(1,i)=newv(1,i);
        v(2,i)=newv(2,i);
        F(1,i)=newF(1,i);
        F(2,i)=newF(2,i);
    end
    leastx=r(1,Ndtot+1);
    for j=0:1:Npb-2
        if leastx<=r(1,Ndtot+2+j)
           leastx=leastx;
        else
           leastx=r(1,Ndtot+2+j);
        end
    end
    mostx=r(1,Ndtot+1);
    for k=0:1:Npb-2
        if mostx>=r(1,Ndtot+2+k)
            mostx=mostx;
        else
            mostx=r(1,Ndtot+2+k);
        end
    end
    extx=mostx-leastx;
    if abs(extx)>LX1
        extx=mostx+leastx;
    end
    fid = fopen('ext.txt', 'a'); % Opening output file
    fprintf(fid,'%-07.4f %-07.4f\r\n',t, extx); %writing value file
    fclose(fid); %Closing output file
    ry=r(2,Ndtot+1:Ntot);
    fid = fopen('tvsrneg13.txt', 'a');
    fprintf(fid,'%-4.2f\r\n',t);
    fprintf(fid,'%-07.4f\r\n',ry);
    fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the particle movement through the channels
    figure(1)
    plot(r(1,N+1:N+Nwall),r(2,N+1:N+Nwall),'o','LineWidth',0.2, ...
    'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b')
    hold on
% whitebg('white')
% set(gcf,'Color',[0.5,1,0.6])
    plot(r(1,1:N),r(2,1:N),'o','MarkerSize',6, ...
    'MarkerEdgeColor','k','MarkerFaceColor','y')
% plot(r(1,200),r(2,200),'>','MarkerSize',6, ...
% 'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(r(1,Ndtot+1:Ntot),r(2,Ndtot+1:Ntot),'-ok','MarkerSize',5, ...
    'MarkerEdgeColor','k','MarkerFaceColor','r')
    text(r(1,Ndtot+1),r(2,Ndtot+1),num2str(Ndtot+1))
    text(r(1,Ndtot+2),r(2,Ndtot+2),num2str(Ndtot+2))
    text(r(1,Ndtot+3),r(2,Ndtot+3),num2str(Ndtot+3))
    text(r(1,Ndtot+4),r(2,Ndtot+4),num2str(Ndtot+4))
    text(r(1,Ndtot+5),r(2,Ndtot+5),num2str(Ndtot+5))
% axis tight
    axis([-LX/2 LX/2 -LY/2 LY/2])
    drawnow
    hold off
% Plotting the averaged velocity in each bin over a set time step
    for i=1:1:N
        if r(2,i)<=bnbtm %&& r(2,i)>=(LY2)
            vnum(1)=vnum(1)+1;
            vavgd(1)=vavgd(1)+v(1,i);
        elseif r(2,i)>=bntp %&& r(2,i)<=(LY1)
            vnum(tbins)=vnum(tbins)+1;
            vavgd(tbins)=vavgd(tbins)+v(1,i);
        end
    end
    for p=1:1:bins
        sect=((p)*nybins)-(LY1)+rcw;
        sect2=((p-1)*nybins)-(LY1)+rcw;
        for i=1:1:N
            if r(2,i)<=sect && r(2,i)>sect2
                vnum(p+1)=vnum(p+1)+1;
                vavgd(p+1)=vavgd(p+1)+v(1,i);
            end
        end
    end
    for tm=1:1:10
        if t==(tf/10)*tm
            vavgd2(:)=vavgd(:)./vnum(:);
            vavgd2(isnan(vavgd2))=0;
            vtot=sum(vavgd2,2);
            vnorm=vtot/tbins;
            vovnrm(1,:)=vavgd2(:)/vnorm;
            figure(4);
            plot(nymat(:),vovnrm(1,:),'o','MarkerSize',6, ...
            'MarkerEdgeColor','k','MarkerFaceColor','k')
            xlabel('BINS')
            ylabel('VELOCITY /AVERAGE VELOCITY')
            axis tight
% axis([-15 15 -1 3])
            set(gca,'XMinorTick','on','YMinorTick','on')
            drawnow
            vnormtxt=[nymat; vovnrm];
% open the file with write permission
            fid = fopen('posnorm1.txt', 'a'); % Opening output file
            fprintf(fid,' \r\n');
            fprintf(fid,'%-4.2f\r\n',t);
            fprintf(fid,'%-07.4f %-07.4f\r\n',vnormtxt); %writing to output file
            fclose(fid); %Closing output file
        else
            continue
        end
    end
    if t<=10
        vnum(1,1:tbins)=0;
        vavgd(1,1:tbins)=0;
    end
end
figure(2);
plot(nymat(:),vavgd2(:),':k','LineWidth',2)
xlabel('BINS')
ylabel('AVERAGE VELOCITY')
axis tight
% axis([-1.5 1.5 -1.5 1.5])
set(gca,'XMinorTick','on','YMinorTick','on')
toc
