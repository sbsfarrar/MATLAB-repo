function [ fd,fs ] = FitMultiGaussScew( xf, Ifit, Lb, Ub, Sp, NG,Nr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i=1:Nr

fo=fitoptions('Method','NonlinearLeastSquares','Lower',Lb,'Upper',Ub,'StartPoint',Sp);

switch(NG)
    
    case(5)
        ft=fittype('A1*exp(-(x-b1)^2/c1^2)*(1+erf(a*(x-b1)/c1))+A2*exp(-(x-b2)^2/c2^2)*(1+erf(a*(x-b2)/c2))+A3*exp(-(x-b3)^2/c3^2)*(1+erf(a*(x-b3)/c3))+A4*exp(-(x-b4)^2/c4^2)*(1+erf(a*(x-b4)/c4))+A5*exp(-(x-b5)^2/c5^2)*(1+erf(a*(x-b5)/c5))','options',fo);
        [fd,fs]=fit(xf,Ifit,ft);
    case(4)
        ft=fittype('A1*exp(-(x-b1)^2/c1^2)*(1+erf(a*(x-b1)/c1))+A2*exp(-(x-b2)^2/c2^2)*(1+erf(a*(x-b2)/c2))+A3*exp(-(x-b3)^2/c3^2)*(1+erf(a*(x-b3)/c3))+A4*exp(-(x-b4)^2/c4^2)*(1+erf(a*(x-b4)/c4))','options',fo);
        [fd,fs]=fit(xf,Ifit,ft);
    case(3)
        ft=fittype('A1*exp(-(x-b1)^2/c1^2)*(1+erf(a*(x-b1)/c1))+A2*exp(-(x-b2)^2/c2^2)*(1+erf(a*(x-b2)/c2))+A3*exp(-(x-b3)^2/c3^2)*(1+erf(a*(x-b3)/c3))','options',fo);
        [fd,fs]=fit(xf,Ifit,ft);
    case(2)
        ft=fittype('A1*exp(-(x-b1)^2/c1^2)*(1+erf(a*(x-b1)/c1))+A2*exp(-(x-b2)^2/c2^2)*(1+erf(a*(x-b2)/c2))','options',fo);
        [fd,fs]=fit(xf,Ifit,ft);
    case(1)
        ft=fittype('A1*exp(-(x-b1)^2/c1^2)*(1+erf(a*(x-b1)/c1))','options',fo);
        [fd,fs]=fit(xf,Ifit,ft);
        
end

Sp=coeffvalues(fd);

end

