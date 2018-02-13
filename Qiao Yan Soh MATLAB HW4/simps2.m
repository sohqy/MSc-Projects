function [int] = simps2(fhandle, lowerlimit, upperlimit )
%SIMPS2 calculates the definite integral of a function between an upper
%and lower limit, using the simpsons rule. 
%
%   FILE NAME: simps2.m
%   The number of intervals is minimized to produce an answer with 4 digit
%   precision through comparing the difference in integration value
%   calculated and that when using twice the number of intervals. 
%
%   INPUTS:     fhandle function, xmin and xmax integration limits
%   OUTPUTS:    Value of definite integral of the function
%   AUTHOR:     Qiao Yan Soh
%   DATE:       1 March 2017

%%
% IF checks for correct input sequence with the upper and lower limits. 
if upperlimit > lowerlimit
    N = 2;                          % Smallest number of intervals we can use to apply Simpsons' Method
    h = (upperlimit-lowerlimit)/N;  % Calculates the width of each interval
    x = lowerlimit:h:upperlimit;    % Generates the list of points at which the fxn should be evaluated at 
    int0 = 1;                       % Initial 'guess' of the integral value
    int = (h/3)*(feval(fhandle, x(1))+feval(fhandle, x(end))+4*feval(fhandle,x(2)));    % Simpsons rule for the least number of intervals.
    
    % WHILE loop is used to minimize the number of intervals required to
    % perform the calculation to the desired precision. 
    while abs(int-int0)/int>=0.0001         % Precision of 4decimal places
        int0 = int;                     % Saves the old calculated integral value into 
        N = 2*N;                        % Increases the number of intervals used
        n = 1:N+1;                      % Creates a list of element index numbers
        h = (upperlimit-lowerlimit)/N;  % Recalculates the width of each interval
        x = lowerlimit:h:upperlimit;    % Creates a new list of evaluation points
        
        coeff = 3 + (-1).^n;            % Creates a list of alternating coefficients 
        coeff(1) = 1;                   % Fixes the coefficient of f1 and fn+1 to 1
        coeff(end) = 1;
    
        f(n) = feval(fhandle,x(n));     % Evaluates the function at each point
        int = (h/3)*sum((coeff.*f));    % Integration calculation through implementing Simpsons' Equation
    end
    
    fprintf('Number of required intervals is ')     % Prints a text header for displaying the interval number
    fprintf('%u\n', N)
else 
    disp('ERROR, wrong sequence input')           % Displays an error when input sequence is incorrect. 

end

