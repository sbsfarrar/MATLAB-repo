function [ Sp,Lb,Ub ] = AssignBoundsScew( delta, NG, TempX, TempI, bt, ct, scew )
%UNTITLED4 Create vectors of startpoints and upper and lower bound  for fit
%   Detailed explanation goes here

Sp=[];
Lb=[];
Ub=[];
AL=[];
bL=[];
cL=[];
aL=[];
AU=[];
bU=[];
cU=[];
A=[];
c=ones(1,NG)*ct;

for i= 1:NG
    index=find(TempX==bt(i));
    A(i)=TempI(index);
    
   
   
   
  
    
end

 %Define Lower bounds
    AL=A-delta(1);
    bL=bt-delta(2);
    cL=c-delta(3);
    %Define Upper bounds
    AU=A+delta(1);
    bU=bt+delta(2);
    cU=c+delta(3);
    
Ub=[AU,scew(2),bU,cU];
Lb=[AL,scew(3),bL,cL];
Sp=[A,scew(1),bt,c];

end

