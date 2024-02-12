function [ Areas, Errors ] = FindScewAreas( FitData, NumPeaks, Norm,Xmax,C )
%Commputes the Areas of the Scewed Guassians and the error in the areas
%based on the confidence interval of the fit parameters.
%FitData is the cell array of fit parameters( usually fd from Fitting
%script)
%NumPeaks is the number of peaks in the profile
%Xmax is the maximum distance to commpute the area ( usually this should be
%large enough so that the right most peak is within Xmax)
%Norm=1 to Normalize the areas, Norm=0 does not normalize the arieas
%C is the value for the confidence interval on the fit parameters
%   Detailed explanation goes here


    for i=1:length(FitData)
        Fnow = FitData{i};
        conf=confint(Fnow,C);
        conf=conf(:,1:NumPeaks);
        for j=1:NumPeaks
            %eval(['Areas(i,j)= Fnow.a', num2str(j) , '*Fnow.c', num2str(j), ';'])
            eval(['A= Fnow.A', num2str(j),';'])
            eval(['b=Fnow.b',num2str(j),';'])
            eval(['c=Fnow.c',num2str(j),';'])
            a=Fnow.a;
            x=0:Xmax;
            cL=conf(1,j);
            cH=conf(2,j);
            w=0.06;
            c1=(1+w)*c;c2=(1-w)*c;
            y=A*exp(-(x-b).^2/c^2).*(1+erf(a*(x-b)/c));
            Atemp=sum(y);
            switch(A<0)
                case(1)
                    %The Area is negative use the absolute value as the
                    %maximum possible, and the minimum to be zero
                    
                    Atemp=abs(Atemp)/2;
                    E=Atemp/2;
                case(0)
                    EH=sum(cH*exp(-(x-b).^2/c1^2).*(1+erf(a*(x-b)/c1)));
                    EL=sum(cL*exp(-(x-b).^2/c2^2).*(1+erf(a*(x-b)/c2)));
                    if EL>0
                        E=(EH-EL)/2;
                    else
%                         EL=Atemp;
%                         E=EL;
                        E=EH/2;
                    end
            end
            

            
            
            Areas(i,j)=Atemp;
            Errors(i,j)=E;
        end
    end
   
    switch (Norm)
        case(1)
            
             for i=1:length(Areas(:,1))
                 T=sum(Areas(i,:));
                 Etot=sum(Errors(i,:))./T;
                Errors(i,:)=sqrt((Errors(i,:)./Areas(i,:)).^2+Etot^2);
                Areas(i,:)=Areas(i,:)./T;
                Errors(i,:)=Errors(i,:).*Areas(i,:);
                
               
             end
        end
    
end

