function [ Areas ] = AreaFromPoint(ap,Ic )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Numpeak=length(ap(:,1))/2;

for i=1:Numpeak

    istart=(i-1)*2+1;
    iend=istart+1;
    
    x1=ap(istart,1);
    x2=ap(iend,1);
    Areas(i)=sum(Ic(x1:x2));
    
end




end

