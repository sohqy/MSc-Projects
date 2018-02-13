%%
%Script M-file: bingo.m
%Description: Simulates a bingo round, generating a winning set of 10
% numbers from a pool of unique numbers ranging from 1-99
% Numbers are not replaced after they are picked.
% 
%Author: Qiao Yan Soh
%Date: 23 Feb 2017

%%

set = 1:99;            % Creates an array of available digits to pick from
win = zeros(1,10);     % Creates an array to store the winning numbers

% FOR loop iterates the picking process over the number of times a winning
% digit is picked
for i = 1:10                         % Sets the loop number
    pick = randi(length(set),1);     % Picks a random index number from the available digits  
    win(i) = set(pick);              % Adds the picked number into the winning array
    set(pick) = 0;                   % Removes the picked number from the available digits
    set=set(set~=0);                 % Updates the array of available digits to disinclude the previously picked digit.
end

fprintf('The winning numbers are...')
fprintf('\n')
fprintf('%u, ', win(1:end-1))        % Prints the 10 results in a line
fprintf('and')
fprintf(' %u\n', win(end))