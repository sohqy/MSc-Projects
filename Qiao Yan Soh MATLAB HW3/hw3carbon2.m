%%
%Script M-file: carbon.m
%Description: Simulates the decay of radioactive carbon by focusing on the
%fate of 100 carbon atoms over an experiment time of 100 minutes. The
%number of remaining atoms is calculated and recorded after each minute. 
% 
%Author: Qiao Yan Soh
%Date: 23 Feb 2017
%%

decay_rate = 0.0338;        % The probability of an atom decaying in each minute (loop)
atoms = ones(1,100);        % Sets up the atomic view, where an undecayed atom = 1 and decayed atom = 0
undecayedcount = [100 zeros(1,100)];    % Stores the number of atoms that remain after each minute, initializing the first value = 100
minute = 1:100;             % Running the experiment for 100 minutes

% <<INITIAL FOR LOOP CODE>>
% Looks at all 100 atoms, decayed or not, over the entire experiment.
% FOR loop runs the decay count each minute, and updates the values
% accordingly. 
% for m = minute
%     undecayed = atoms.*rand(1,100) > decay_rate;   % rand gives a random number for determining decay, and is set to 0 if the atom has already decayed. 
%     atoms = undecayed;                             % Updates the atomic view of the 100 atoms of interest 
%     undecayedcount(m+1) = sum(undecayed);          % Counts the total number of remaining atoms after this 'minute' and stores it in the counter
% end

% << UPDATED FOR LOOP CODE>>
% Comparatively more accurate - since it stops keeping track of any decayed atoms and
% reduces time by reducing the number of elements it needs to look at. This
% is suitable since the atoms are not distinguishable from each other and
% its element index is unimportant. 
% FOR loop runs the decay count each minute, and updates the values
% accordingly. An atom decays, and value is set to 0 if its corresponding randomly generated
% number is less than the decay rate. 
for m = minute
    undecayed = rand(1, length(atoms)) > decay_rate;    % Produces an array of random numbers, and compares it to the decay rate to set 'decayed' atoms to 0
    atoms = undecayed(undecayed~=0);                    % Updates the atomic view to keep track of remaining atoms only. 
    undecayedcount(m+1) = sum(atoms);                   % Counts the total number of remaining atoms after this 'minute' and stores it in the counter
end

close all
figure(1)
bar(0:100, undecayedcount, 'b')   % Plots the numbers in a bar chart
hold on
axis([0 100 0 100])
xlabel('Time (Minutes)')
ylabel('Number of remaining carbon atoms')
title('Radioactive Carbon Decay Count')
fplot(@(t) 100*exp(-decay_rate*t),[0 100], 'r', 'linewidth', 2)   % Plots the theoretical decay prediction expression
legend('Simulation numbers', 'Theoretical expression')
hold off
