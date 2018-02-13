function [ysto] = mysqrtf(r)
%MYSQRTF calculates the squareroot of a non-negative real number iteratively using
%the Babylonian method.

%   File name: mysqrtf.m
%   The result is a vector array of the calculated squareroot values of the
%   input. 
%   Result is calculated through an initial estimate of the squareroot
%   value, and calculating a second estimate such that their average is a
%   better approximation of the result. This process is iterated until a
%   suitable accuracy (< 10*eps) is achieved. 

%   INPUTS:     Vector array or single value r
%   OUTPUTS:    Vector array or single value ysto
%   Author:     Qiao Yan Soh
%   Date:       21 February 2017

%%
ysto = zeros(1, length(r));     % Creates a zero array of the required length for updating to increase speed, rather than appending a new element each time

if isreal(r)==1                 % Condition for real numbers in given array
    if sum(r < 0) == 0          % Condition for positive numbers in given array
        for i = 1:length(r)     % Labels each element in the array with an index
            for x = r(i)        % Picks out the corresponding element value using the defined index
                 
                if x == 0       % Babylonian method cannot be used to calculate squareroot of 0 
                   ysto(i)=0;   % Sets the result to 0

                else            % Runs the Babylonian method calculations if number is real and non negative
                    y_o = 1;                    % Defines the first initial guess
                    y = (y_o+x/y_o)/2;          % Sets an initial calculation for y, the solution according to the first initial guess

                    % WHILE loop iterates the squareroot approximation until a suitable
                    % accuracy is reached.
                    while (abs(y-y_o)/y > 10*eps)       % Compares the calculated value to the initial guess, to a suitable difference
                        y_o = y;                        % Stores the previously calculated y value into y_o
                        y = (y+x/y)/2;                  % Calculates and stores a new y value based on the previously calculated y value
                    end

                    ysto(i) = y;                % Stores the calculated y value into an array ysto
                end 
            end
        end
        
    % If negative numbers are present, display an error message 
    else
        disp('ERROR: A negative number was provided.')  % Displays an error message
        ysto = 'Array given is in incorrect form';      % Rewrites the function output to show a more obvious error message 
    end
    
% If imaginary numbers are present, display error message 
else
    disp('ERROR: Array contains imaginary numbers')
    ysto = 'Array given is in incorrect form';
end