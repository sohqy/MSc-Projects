%% 
%Script M-file: Erathosthenes3.m
%Description:
% Determines all prime numbers between 1 and a user given n,
% using the sieve of Eratosthenes, printing the result alongside an output
% message. 
%Author: Qiao Yan Soh
%Date: 

%% 

prompt = 'Please provide an integer n: ';   % Asks the user for an input value, storing question as 'prompt'
n = input(prompt);      % Stores the user provided value into n

x =[0 2:n];             % Initializes the set of numbers to test

% Loop sieves out prime numbers from vector x
for k = 2:10            % Sets range of values of k
    for m = 2:100/k     % Creates the required list of multiples
        x(m*k)=0;       % Sets all multiples of k to 0. 
    end
end

y=nonzeros(x);          % Compiles all prime numbers into a single vector

% Display results in text form
fprintf('List of all prime numbers between 1 and ')
fprintf('%u', n)
fprintf(' are: ')
fprintf('%u,',y(1:end-1))
fprintf('%u. \n', y(end))