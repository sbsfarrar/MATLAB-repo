h=10.^(-20:0);
errors = abs(cos(2)-(sin(2+h)-sin(2))./h);
loglog(h,errors,'.-');
grid