   
        
        for j=1:NG
            eval(['A=fdt.A', num2str(j),';']);
            eval(['b=fdt.b',num2str(j),';']);
            eval(['c=fdt.c',num2str(j),';']);
            a=fdt.a;
            x=0:Xmax;
            y=A*exp(-(x-b).^2/c^2).*(1+erf(a*(x-b)/c));
            hold on;
            plot(x,y,'b');
        end

    hold off
    
    
   