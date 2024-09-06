function [b, c, delta,scew] = ChangeBoundsScew(bt,ct,deltat,scewt);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Q=input('Change b(b), c(c), delta (d), scew (s) or all (all): ','s')

switch(Q)
    
        
    case('b')
        s=input('New Vector (n) or delta (d) ?: ','s');
        switch(s)
            case('n')
         b=input('Enter a vector of starting positions for peaks centers [b1 b2...]: ');
        c=ct
        delta=deltat;
        scew=scewt;
            case('d')
                d=input('Enter shift distance: ');
                b=bt+d;
                c=ct
        delta=deltat; scew=scewt;
        end
        
    case('c')
         b=bt;
        c=input('Enter a starting guess for the width: ');
        delta=deltat; scew=scewt;
        
    case('d')
         b=bt;
        c=ct; scew=scewt;
        delta=input('Change in fit parameters absolute change for db [da db dc]: ');
        
    case('all')
        b=input('Enter a vector of starting positions for peaks centers [b1 b2...]: ');
        c=input('Enter a starting guess for the width: ');
        delta=input('Change in fit parameters absolute change for db [da db dc]: ');
        scew=input('Enter a starting guess for skew, Upper, Lower: ');
    case('s')
        scew=input('Enter a starting guess for skew, Upper, Lower: ');
        b=bt;
        c=ct;
        delta=deltat;
        
end

end

