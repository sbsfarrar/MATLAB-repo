function [ Areas ] = ManualAreas( Iall,Xall,Numpeaks )
%ManualAreas 
%   Detailed explanation goes here

ap=[];
Nmax=Numpeaks*2;
for i =1:length(Xall)
    
    Xc=Xall{i};
    Ic=Iall{i};
    figure(1);
    hold off;
    plot(Xc,Ic);
    
    if i==1
        
        
        
        L=length(Xc)
        Yall=0;
        coords(1,1:2)=[1,Yall];
        d=Nmax-1;
        for j=2:Nmax
            coords=[coords;[L*(j-1)/d,Yall]];
        end
        
        bgbox=impoly(gca,[coords],'Closed',false);
        
        input('Draw points for calculating areas ');
        
        BGcoord=round(getPosition(bgbox));
        ap(:,1)=BGcoord(:,1);
        ap(:,2)=Ic(BGcoord(:,1))';
        hold on;
        plot(ap(:,1),ap(:,2),'ro');
        Areas(i,:) = AreaFromPoint(ap,Ic );
    else
        hold on;
        plot(ap(:,1),ap(:,2),'ro');
        check=input('OK? Good(1) Redraw Points(0) ');
        
        if check==1
            Areas(i,:) = AreaFromPoint(ap,Ic);
        else
            input('Draw points for calculating areas ');
            
            BGcoord=round(getPosition(bgbox));
            ap(:,1)=BGcoord(:,1);
            ap(:,2)=Ic(BGcoord(:,1))';
            hold on;
            plot(ap(:,1),ap(:,2),'ro');
            Areas(i,:) = AreaFromPoint(ap,Ic );
        end
        
    end
    
end

