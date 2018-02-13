function [x] = bubblesortf(x)
%BUBBLESORTF sorts an input vector in ascending order

%   File name: bubblesortf.m
%   The result is a sorted array where elements are placed in ascending
%   order.
%   The result is produced by comparing the first two elements, and
%   swapping their positions if the second element is smaller than the
%   first. This is comparison is also performed for incrementing pairs through the vector, 
%   and the entire process is iterated through n times, where n is the length of
%   the vector.

%   INPUTS:    Vector x
%   OUTPUTS:   Vector x with elements in ascending order
%   Author:    Qiao Yan Soh
%   Date:      21 February 2017

%%
% FOR loop changes the stored x array by swapping the positions of each
% element. 
for m = 1:length(x)         % Repeats the process n times
    for n = 1:length(x)-1   % Sequentially compares the values in each pair of elements
        if x(n+1) < x(n)
            x([n n+1]) = x([n+1 n]);    % Swapping element positions
        end
    end
end

end

