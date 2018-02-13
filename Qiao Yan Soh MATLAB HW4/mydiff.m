function [ d1 ] = mydiff(fhandle, x)
%MYDIFF Numerically calculates the first derivative of the input function,
%at the user defined x. 
%   
%   FILE NAME: mydiff.m
%   Routine starts with a step size of 1, and increments by a order of
%   magnitude each time the calculated differential value does not fulfil
%   the precision requirement. 

%   INPUTS:     fhandle function, single value x
%   OUTPUTS:    Calculated value of the first derivative
%   AUTHOR:     Qiao Yan Soh
%   DATE:       1 March 2017

%%
h = 1;      % Initial step size
d0 = 1;     % Initial guess of the differential value
d1 = (feval(fhandle, x+h)-feval(fhandle, x))/h;     % First principles formula for differentiation

 % WHILE loop is used to optimize the step size required to perform the
 % calculation to a desired precision. 
while abs(d1-d0)/d1 >= 0.0001       % Precision of 4sf required
    d0 = d1;                        % Updating the old differential value for future comparison
    h = h/10;                       % Decreasing the step size
    d1 = (feval(fhandle, x+h)-feval(fhandle, x))/h;     % Recalculating the differential according to the new step size
end

fprintf('Required step size: ')
fprintf('%u \n', h)


