%%

%% Calculating Unperturbed Wavefunctions 
a = 100;                        % Length of square well
lowerlimit = -a/2;              % Center of well at position x = 0
upperlimit = a/2;               % Calculating boundary positions
x = linspace(-80, 80, 1000);    % Create an array of positional values

% Normalisation constant is the same for all functions in an infinite
% square well. 
phi1 = @(x) cos((pi*x)/a).^2;   % Using the ground state for calculating the normalisation constant
norm = 1/sqrt(integral(phi1, lowerlimit, upperlimit))   % Normalisation constant

V0 = rectangularPulse(-50, 50,x);       % Defining the infinite square well (Implementing BCs)
fxn_unpert = zeros(10, length(x));      % Generates 10 rows, each to calculate wavefxn values at each x position

for N = 1:10            % Calculating 10 unperturbed wavefunctions, using given forms
    if mod(N,2)==0      % Checking for even indices
        fxn_unpert(N,:) = V0.*norm.*sin((N.*pi.*x)/a);  % Applying normalisation constant and BCs
    else
        fxn_unpert(N,:) = V0.*norm.*cos((N.*pi.*x)/a);
    end
end

% Plotting the first 5 wavefunctions, distinguishing between odd and even
% levels. 
figure(1)
plot(x, fxn_unpert(1:2:5,:), x, fxn_unpert(2:2:5,:),'--')
legend('1','3','5','2','4')
title('Visualisation of the first 5 eigenfunctions in an infinite square well')
ylabel('Normalised Eigenfunction Amplitude')
xlabel('Position x (nm)')
axis([-80,80, -0.2, 0.2])

%% Adding Pertubation to Well
delta = 5*10^-22*10^18;                     % Pertubation in units of nm
V1 = rectangularPulse(a/4, a/2, x)*delta;   % Defining the pertubation potential

hbar = 1.05*10^-34*10^18;   % In units of nm^2 kg s^-1
m_e = 9.1*10^-31;           % Units of kg

% Calculating the unperturbed energies
E_unperturb = zeros(1,10);      % Initialising array for energies
for n = 1:10                    % Calculates unperturbed energies for each state
    E_unperturb(n) = (pi*hbar*n)^2/(2*m_e*a^2);
end

H0 = diag(E_unperturb);     % By definition of unperturbed hamiltonian matrix

% Finding Matrix form of pertubation
for n = 1:10
    for m = 1:10
        H1(n,m) = trapz(x, V1.*fxn_unpert(n,:).*fxn_unpert(m,:));
    end
end 

% Use inbuilt eig function to solve perturbed system
H_total = H0 + H1;          % Total hamiltonian of the perturbed system
[eigenfxn, eigenvalues] = eig(H_total);     % Finding eigenfunctions and eigenenergies of the perturbed system
fxn_perturb = -eigenfxn'*fxn_unpert;        % Superpositioning the wavefunction calculations

E_perturb = max(eigenvalues);               % Flatten diagonal matrix of eigenenergies
GroundE_i = find(E_perturb == min(E_perturb));  % Find index of minimum energy
E_perturb(GroundE_i) = 1;                       % Make index value large to disinclude from next find function
FirstE_i = find(E_perturb == min(E_perturb));   % Find index of second smallest energy state

% Plot the two lowest energy states for the perturbed and unperturbed
% system to the same axis for comparison. 
figure(2)
plot(x, fxn_unpert(1:2,:))
hold on
plot(x, fxn_perturb(FirstE_i,:), x, fxn_perturb(GroundE_i,:))
legend('Unperturbed ground state', 'Unperturbed first excited state', 'Perturbed ground state', 'Perturbed first excited state', 'Location', 'southeast')
line([a/4 a/4], [-0.2 0.2],'LineStyle', ':', 'Color','k')       % Drawing perturbation boundary
line([a/2 a/2], [-0.2 0.2],'LineStyle', '--', 'Color','k')      % Drawing well boundaries
line([-a/2 -a/2], [-0.2 0.2],'LineStyle', '--', 'Color','k')
hold off
title('Comparison between the two lowest energy wavefunctions')
ylabel('Normalised Eigenfunction Amplitude')
xlabel('Position x (nm)')
axis([-80,80, -0.2, 0.2])