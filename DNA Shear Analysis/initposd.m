function [ r ] = initposd( N,Nwall )
    global LX LY Np Npb Ntot
    global Q Qm1 nx ny nwallx Nwall2
    global LX2 LY2 lseg nyp
%Initiate Arrays
    ri=zeros(2,Q); r=zeros(2,Ntot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluid Particle Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    for s=0:2:Qm1-1
        ri(1,s+1)=(s*nx)-(LX/2)+0.05;
        if ri(1,s+1)>(LX/2)
            break
        end
        for p=0:2:Qm1-1
            ri(2,p+1)=(p*ny)-(LY/2)+0.05;
            if ri(2,p+1)>(LY/2)
                break
            end
        i=i+1;
        r(1,i)=ri(1,s+1); % X-Position of the particles
        r(2,i)=ri(2,p+1); % Y-Position of the particles
        end
    end
    for s=1:2:Qm1
        ri(1,s+1)=(s*nx)-(LX/2)+0.05+(nx/20);
        if ri(1,s+1)>(LX/2)
            break
        end
        for p=1:2:Qm1
            ri(2,p+1)=(p*ny)-(LY/2)+0.05+(ny/8);
            if ri(2,p+1)>(LY/2)
                break
            end
            i=i+1;
            r(1,i)=ri(1,s+1); % X-Position of the particles
            r(2,i)=ri(2,p+1); % Y-Position of the particles
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wall Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=N+1:1:(N+Nwall2)
        r(1,i)= (i*nwallx)-(LX/2)+0.001-((N+1)*nwallx);
        if r(1,i)>(LX/2)
            break
        end
        r(2,i)=-(LY/2);
        if r(2,i)>(LY/2)
            break
        end
    end
    for i=N+1:1:(N+Nwall2)
        r(1,Nwall2+i)= (i*nwallx)-(LX/2)+0.001-((N+1)*nwallx);
        if r(1,Nwall2+i)>(LX/2)
            break
        end
        r(2,Nwall2+i)=(LY/2);
        if r(2,Nwall2+i)>(LY/2)
            break
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNA Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=N+Nwall;
    rpi(1,1)=LX2;
    for p=1:1:Np
        rpi(2,p+1)=LY2+((p-1)*nyp)+(nyp/2);
        for s=1:1:Npb
            rpi(1,s+1)=rpi(1,s)+lseg;
            if rpi(1,s+1)>(LX/2)
                break
            end
            i=i+1;
            r(1,i)=rpi(1,s+1); % X-Position of the particles
            r(2,i)=rpi(2,p+1); % Y-Position of the particles
        end
    end
end