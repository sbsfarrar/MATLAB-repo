function [ xf,xnf, Ifit, Infit ] = DefineFitRegion( fr,TempI )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



Nr=length(fr)/2;
xf=[];
xnf=[];
xnfb=1;
for i=1:Nr
    j=2*i-1;
    xs=fr(j);
    xe=fr(j+1);
    xnf=[xnf,xnfb:xs-1];
    xf=[xf,xs:xe];
    xnfb=xe+1;
end

Ifit=TempI(xf);
Infit=TempI(xnf);

end

