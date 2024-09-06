% C=0;
% Mod=0;
% Iall=TempIall;
% Xall=TempXall;
% Mod=0;
for i=Start:length(Iall)
%% 

%% Load Data

TempI=Iall{i,1};
TempX=Xall{i};
hold off
plot(TempX,TempI);
%% 
C=0;
while C==0
%% Input Fit Parameters
if (i==1 && Mod==0)
%         Pn=input('Enter peak number to use as reference: ');
%         Rbound=input('Enter Right distance from peak to integrate: ');
%         Lbound=input('Enter Left distance from peak to integrate: ');
        
        Pn=1;
        Rbound=1;
        Lbound=1;
        
        
        NG=input('Enter number of peaks to fit with Gaussian: ');
        display('Format is a*exp(-((x-b)/c)^2)*ScewFunction');

        bt=input('Enter a vector of starting positions for peaks centers [b1 b2...]: ');
        ct=input('Enter a starting guess for the width: ');
        scew=input('Enter a starting guess for skew, Upper, Lower: ');
        delta=input('Change in fit parameters absolute change for [da db dc]: ');
        fr=input('Enter vector desribing fit region in even pairs [xb1 xe1 xb2 xe2..]: ');
        Nr=input('Enter number of refits: ');
else
    if Mod>0
        Q=input('Change bounds (b) Change region (r) Change Number Peaks(p): ','s');
        plot(TempX,TempI);
        switch(Q)
            case('b')
                [bt, ct, delta,scew]=ChangeBoundsScew(bt,ct,delta,scew);
            case('r')
                fr=input('Enter vector desribing fit region in even pairs [xb1 xe1 xb2 xe2..]: ');
            case('p')
                 NG=input('Enter number of peaks to fit with Gaussian: ');
        end
    end
end
%% Assign Bounds
Sp=[];
Lb=[];
Ub=[];

[ Sp,Lb,Ub ] = AssignBoundsScew( delta, NG, TempX, TempI, bt, ct, scew );
%% %define fit region



[ xf,xnf, Ifit, Infit ] = DefineFitRegion( fr,TempI );

%% Plot Defined Region
hold off;
plot(xnf,Infit,'ro');hold on;plot(xf,Ifit,'bo');
%% Fit Data

[ fdt,fst ] = FitMultiGaussScew( xf', Ifit', Lb, Ub, Sp, NG,Nr );
hold on;

plot(fdt);
Af=[];
bf=[];
for k=1:NG
    eval(['Af(k)=fdt.A',num2str(k),';']);
     eval(['bf(k)=fdt.b',num2str(k),';']);
    
    
end
plot(bf,Af,'go','MarkerFaceColor','g');
Xmax=max(xf);
SingleGaussPlotFromMultiFitScew
CI=confint(fdt,0.95);
 f=find(CI(1,:)~=CI(1,:));
 if (isempty(f)==0)
     display('Fit Parameters Fixed at Bound');
     
 end
 fdt
C=input('Is fit acceptable? (1 or 0): ');
p=mean(CI);
p=p(2:3:end-1);

if C==0
    Mod=1;
else 
    Mod=0;
%     bt
%     eval(['pL=fdt.b',num2str(Pn),';']);
%     pL=round(pL);
%     Istart=pL-Lbound;
%     Iend=pL+Rbound;
    %Itot(i)=sum(TempI(Istart:Iend));
    
    fd{i}=fdt;
end

end
%% 

end

%% 




