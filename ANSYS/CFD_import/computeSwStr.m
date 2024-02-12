function lambda_ci = computeSwStr(dxux, dxuy, dxuz, ...
                                  dyux, dyuy, dyuz, ...
                                  dzux, dzuy, dzuz)
                              
%% Copyright 2015 The MathWorks, Inc.

% This function is called by CFD_import.m and computes the swirling
% strength (vorticity) based on the Lambda Ci cirterion.

%%                              
% build the velocity tensor
velocity_tensor = [dxux, dxuy, dxuz;
                   dyux, dyuy, dyuz;
                   dzux, dzuy, dzuz];
               
% compute the eigenvalues               
eigenvalues = eig(velocity_tensor);

% is eigvalues are complex set the value
if ~isreal(eigenvalues)
    lambda_ci = sum(abs(eigenvalues))/2;
else
    lambda_ci = 0;
end
