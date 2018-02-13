%% SCRIPT INFORMATION
%
% Script fxn m-file:    TDSE2.m
% Description:  
%   Solves the Time Dependent Schrodinger equation to determine
%   the time evolution and simulate an electron moving through a potential,
%   through the Finite Domain Time Difference method. 
%   Produces an animation of the wavefunction and probability density
%   evolution. 
%
% Author:   Qiao Yan Soh (00866106)
% Date:     22 March 2017

close all

%% Standard values
m_e = 9.109*10^-31;             % Mass of electron in kg
hbar = 1.05*10^-34*10^18;       % Planck's Constant in nm^2 kg s^-1
h = hbar*2*pi;
e = 1.9e-16;                    % Elementary charge

%% Set up of Simulation
% --- Simulation canvas options:
L = 500;                        % Simulation window, in nm
n_gridpts = 1000;               % Number of simulation spatial divisions
n_timestep = 15000;             % Number of simulation time divisions
framestep = 100;                % Output animation jump frames, in seconds

% Internal set up calculations of canvas
x = linspace(0,L, n_gridpts);   % Setting up x-axis 
dx = x(2)-x(1);                 % grid spatial step, in nm

% --- Physical picture options:
% Potential options
% Default options:  wellheight = 100eV; wellwidth = 200nm; 
%                   well center = simulation window centre = n_gridpts/2
wellheight = 80;                   % Amplitude/Height of well in eV
wellwidth = 200;                    % Width of well in nm
potential = 'quadratic';           % Choose well type: 'free', 'square','quadratic'
well_centre_index = n_gridpts/2;    % The element index of the centre of the simulation canvas

% Particle wave options
% Default options:  wavelength = 30nm; wavepacket width = 20nm
wavelength = 30;                % in nm
s = 20;                         % Wavepacket width
wavecentre = 250;               % User input of wavecentre location
wavecentreindex = find(round(x) == wavecentre, 1);      % x is not in integer increments, thus requiring round. 
k = 2*pi/wavelength;            % Wavevector

%% Potential well set up
V = zeros(1, length(x));        % Initialize potential

% Switch carries out the code under the specified case when the 
% switch_expression is set to the case_expression. 
switch potential
    case 'free'             % No potential, free particle
        
    case 'quadratic'        % Quadratic/Harmonic trap
        V1 = rectangularPulse(L/2-wellwidth/2, L/2+wellwidth/2, x); % Imposing well boundaries
        quadfxn = @(x) (x-L/2).^2/10000;                            % Creating the basic well shape
        shift = feval(quadfxn, L/2-wellwidth/2);                    % Calculating the amount of shifting required for a well
        V = wellheight*V1.*(quadfxn(x)-shift);                      % Adjusting the well according to the required well height and boundaries
        
    case 'square'           % Square well
        lower_wall_index = well_centre_index - round(wellwidth/2);  % Determining corresponding element
        upper_wall_index = lower_wall_index + wellwidth;
        V(lower_wall_index:upper_wall_index) =  - wellheight;       % Setting the range of elements to the correct potential
    
    case 'doublequad'       % Double quadratic wells - for determining tunnelling possibilities
        V1 = rectangularPulse(L/2-wellwidth/2, L/2+wellwidth/2, x);         % Setting well boundaries
        quadfxn = @(x) ((x-L/2).^4-(wellwidth/2*(x-L/2)).^2)/(2.5*10^7);    % Creating the basic well shape
        V = wellheight.*V1.*quadfxn(x);                                     % Adjusting the final potential well according to the required height

    case 'step'           % Square well
        lower_wall_index = well_centre_index - round(wellwidth/2);  % Determining corresponding element
        upper_wall_index = lower_wall_index + wellwidth;
        V(lower_wall_index:upper_wall_index) =   wellheight;        % Setting the range of elements to the correct potential
end

%% Initial wavefxn
phi_r = zeros(1, length(x));        % Initializing wave array that stores 
phi_i = zeros(1, length(x));        % wave values at each gridpoint

% Calculate the initial value of each gridpoint. 
for nx = 1:n_gridpts
    phi_r(nx) = exp(-0.5*((x(nx)-x(wavecentreindex))/s).^2)*cos(k*(x(nx)-x(wavecentreindex)));
    phi_i(nx) = exp(-0.5*((x(nx)-x(wavecentreindex))/s).^2)*sin(k*(x(nx)-x(wavecentreindex)));
end
% Normalizing wavefunction
norm = trapz(phi_r.^2 + phi_i.^2);  % Determining the normalisation constant of the wavepacket
phi_r = phi_r./sqrt(norm);          
phi_i = phi_i./sqrt(norm);

% Saving initial wavefunction values
phi_i1 = phi_i;
phi_r1 = phi_r;

%% Visualisation of initial set up
% Figure 1 is static- displays and compares wavefunction evolution.
% This set of code produces the initial set up of the system. 
figure(1)                       
subplot(2,1,1)                  
yyaxis right
plot(x,V)                       % Plot double y-axes to visualise well positions. 
xlabel('Position x (nm)')
ylabel('Potential (eV)')

yyaxis left
plot(x, phi_i1, 'r', x, phi_r1, 'b')    % Plots initial wavefunctions
title('Initial set up of simulation')
ylabel('Normalised \phi')
legend('\phi_{im}', '\phi_{re}')

% Figure 2 is animated, but has an overlay of the potential function for
% better visualisation.
% This set of code sets up the potential overlay.
figure(2)                       
subplot(2,1,1)                  
yyaxis right
plot(x,V, ':')                   
axis([0, L, -100, 0])
ylabel('Potential (eV)')
xlabel('Position x (nm)')

subplot(2,1,2)
yyaxis right
plot(x,V, ':')
ylabel('Potential (eV)')
axis([0, L, -100, 0])

%% FDTD Propagation of wave

% --- Required FDTD Constants
c1 = 1/10;              % Determines propagation speed of free particle
dt = c1 * 2 * m_e * (dx)^2 * 10^11 / hbar;  % Calculated time step
c2 = e  * dt/ hbar;     % Determines propagation speed of particle in a potential

% --- FDTD Calculation
% FOR loop calculates and appends the wavefunction after each time
% increment step. 
for nt = 1:n_timestep
    % Internal FOR loops calculates the value of each gridpoint
    % individually.
    for nx = 2:(n_gridpts - 1)
        phi_r(nx) = phi_r(nx) - c1*(phi_i(nx+1)-2*phi_i(nx)+phi_i(nx-1)) + c2*V(nx)*phi_i(nx);
    end
    
    for nx = 2:(n_gridpts -1)
        phi_i(nx) = phi_i(nx) + c1*(phi_r(nx+1)-2*phi_r(nx)+phi_r(nx-1)) - c2*V(nx)*phi_r(nx);
    end
    
    Probability = phi_r.^2 + phi_i.^2;      % For each slice
    Timelapsed = ['Time Lapsed = ' num2str(nt*dt, '%2.2f') ' s'];
    fcount = rem(nt, framestep);
    % --- Animation Production
    if rem(nt, framestep)== 0       % To speed up animation/simulation process by jumping frames
        figure(2)                   % Animation figure, to show evolution.
        subplot(2,1,1)
        yyaxis left
        plot(x, phi_i,'r', x, phi_r, 'b')
        title('Wavefunction Evolution')
        axis([0, L, -0.2,0.2])
        text(350, 0.15, Timelapsed)
        ylabel('\phi')
        
        subplot(2,1,2)
        yyaxis left
        plot(x, Probability)
        title('Probability Evolution')
        xlabel('Position x (nm)')
        ylabel('Probability')
        axis([0,L, 0, 0.05])   
        
        pause(0.001)        % Displays frame for a period of time. 
        set(gcf,'color','w');
        drawnow;            % Redraws and updates the frame
        frame = getframe(2);
        im = frame2im(frame);
        [imind, cmap] = rgb2ind(im, 256); 
        outfile = 'TDSE_trapquad.gif';
        
        if nt==100
            imwrite(imind,cmap,outfile,'gif','DelayTime',0, 'Loopcount', inf);
        else
            imwrite(imind,cmap,outfile,'gif','DelayTime',0,'WriteMode','Append');
        end
    
    end
   
end

% Updates the figure to show a static comparison between the initial and
% evolved system. 
figure(1)
subplot(2,1,2)
yyaxis right
plot(x,V)
xlabel('Position x (nm)')
ylabel('Potential (eV)')

yyaxis left
plot(x, phi_i, 'r', x, phi_r, 'b')
title('Evolved simulation')
ylabel('Normalised \phi')
legend('\phi_{im}', '\phi_{re}')    
    
