%% Problem 1: Erathostenes3
erathostenes3

%% Problem 2: Bubblesort

v = randi(100,[1,20]);      % Generates an 20 element array of random integers
y = v;                      % Stores the array in its original sequence in y
indices = 1:length(v);      % Creates an array of index numbers
bsort = bubblesortf(v);     % Calls the written function file, which edits v
      
y1 = sort(y);               % Calls the built-in sort function

figure(1)                   % Visual comparison of the two sorting functions
plot(indices, y, 'ko', indices, bsort, 'rx', indices, y1,'b')   % Plots everything on the same graph
title('Visualisation for Problem 2')
legend('Initial Sequence', 'Written Programme','Built-in Programme','location','best')
xlabel('Index Number')
ylabel('Element Value')

%% Problem 3: MySqrt

a = linspace(0,10);         % Generate 100 linear points between 0 and 10

test = mysqrtf(a);          % Calls the written function file, saving the results into test
bi = sqrt(a);               % Calls the built-in sqrt function, saving the results into bi

figure(2)                   % Visual comparison between the two functions
plot(a, test, 'o', a, bi, 'x')
title('Visualisation for Problem 3')
legend('Written Programme', 'Built-in Programme','location', 'southeast')
xlabel('x')
ylabel('$\sqrt{x}$', 'Interpreter', 'latex')