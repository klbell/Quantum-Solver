function [I, T] = opvlevels()
    clf;
    %% User inputs
    % Energies in the region, widths of regions in nm
    levels = [-4.7, 100; ...
              -5.2, 20; ...
              -4.3, 50; ...
              -3.31, 50; ...
              -3, 0.8; ...
              -2.7, 100];
    
    % Applied external voltage
    voltage = 0.8;
    
    % bandGap in eV
    bandGap = 1;
    
    % Quality
    N = 1000; % Positive interger, number of spatial divisions
    
    % Total time to run for, in femtoseconds
    simulationTime = 10000;
    
    % Timestep, in femtoseconds
    dt = 5;
    
    %% Physical constants
    nm = 1e-9; % meters
    hbar = 4.1356e-15; % eV*s
    me = 9.11e-31; % kg
    q = 1.602e-19; % elementary charge
    
     %% Read data for AM1.5G spectrum
    wavelength = open('wavelength.mat');
    irradiance = open('irradiance.mat');
    wavelength = wavelength.wavelength;
    irradiance = irradiance.irradiance;
   
    %% Potential and initial wavefunction
    % Potential energy generated from user inputs
    [U, structureSize] = generatePotential(levels, voltage, N);
    
    % x, k , t axis
    dt = dt * 1e-15;
    simulationTime = simulationTime * 1e-15;
    dk = 2*pi/structureSize; % Momentum step
    km = N * dk; % Max momentum
    ka = 20; % Momentum Spread
    x = linspace(0, structureSize, N)';
    k = linspace(0, +km, N)';
    
    % Spatial discretness
    dx = x(2)-x(1);
    
    % Initial wavefunctions
    electrons = struct('ID', [], 'alive', [], 'psi', [], 'evolve', []);
    electrons = generateElectron(electrons, randi([floor(0.4*N),floor(0.6*N)],1,1), 2, structureSize, N, U, dt, hbar);
    electrons.alive(1) = 0;
    
    %% Simulate
    % Loop through particles and operate with time evolution operator
    iterations = simulationTime/dt;
    
    % Current
    IL = 0;
    IR = 0;
    I = 0;
    for t = 1:iterations
        
        % How many particles exist
        nParts = countParts(electrons);
        
        % Evolve each particle and add the wavefunctions
        rho = zeros(N, 1);
        for n = 1:length(electrons.ID)
            
            % So far, 5% of the time, generate a new electron
            if (randi(1000,1,1) >= 1000)
                newEnergy = getRandomEnergyPhoton(wavelength, irradiance, hbar);
               
                % Check bandgap
                if newEnergy > bandGap
                    newEnergy = newEnergy - bandGap;
                    electrons = generateElectron(electrons, randi([floor(0.4*N),floor(0.6*N)],1,1), newEnergy, structureSize, N, U, dt, hbar);
                end
                
            end
            if electrons.alive(n) == 1
                % Individual wavefunction evolution
                electrons.psi(:,n) = electrons.evolve(:,:,n)*electrons.psi(:,n);

                % Add the individual probability density
                rho = rho + electrons.psi(:,n).*conj(electrons.psi(:,n));

                % Check if the particle has made it to a contact
                [L, R, destroy] = checkContacts(electrons.psi(:,n), levels, structureSize, N);
                
                % If reached the contacts, kill the particle
                if(destroy)
                   electrons.alive(n) = 0; 
                end
                
                % Calculate left and right electron flow
                IL = IL + q*L;
                IR = IR + q*R;
                I = (IL - IR) / t / dt;
            end
        end
        
        clc; 
        fprintf('Total current = %12.7f microamps \nElapsed Time = %12f ps', I/(1e-6), t*dt / (1e-12));
        % Plot result on a superimposed potential
        %multiPlot(1, x, rho, x, U, 0.1);
    end
    T = iterations*dt;
end
%% Check contacts
% When a particle reaches the contacts, kill it
function [L, R, destroy] = checkContacts(wavefunction, energy, structureSize, N)

    destroy = false;
    L = false;
    R = false;
    leftThreshold = N * energy(1, 2) / structureSize / 2;
    rightThreshold = N * (structureSize - energy(length(energy), 2) / 2)/structureSize;
    
    [~, I] = max(wavefunction);
    
    if (I < leftThreshold)
        destroy = true;
        L = true;
    elseif ( I > rightThreshold)
        destroy = true;
        R = true;
    end
end

%% Count particles
% Count the number of particles that still exist
function number = countParts(electrons)
    possibleNum = length(electrons.alive);
    number = 0;
    for i = 1:possibleNum
        if electrons.alive(i) == 1
            number = number + 1;
        end
    end
end

%% Generate Electron
% Creates a new electron at a given location and energy
function electrons = generateElectron(electrons, location, energy, structureSize, N, U, dt, hbar)
    dk = 2*pi/structureSize; % Momentum step
    km = N * dk; % Max momentum
    x = linspace(0, structureSize, N)';
    dx = x(2) - x(1);
    k = linspace(0, +km, N)';
    k0 = sqrt(2 * energy);
    ka = 10*k0; % Momentum Spread
    const = location/2/k0;

    % Check if there is a free spot
    available = false;
    newLoc = 0;
    for i = 1:length(electrons.ID)
        if electrons.alive(i) == 0
            available = true;
            newLoc = i;
        end
    end
    
    % Get the new ID
    if available
        newID = newLoc;
    else
        numElectrons = length(electrons.ID);
        newID = numElectrons + 1;
    end
    
    % Append the new ID
    if available
        electrons.alive(newID) = 1;
    else
        electrons.ID = [electrons.ID, newID];
        electrons.alive(electrons.ID(newID)) = 1;
    end
    
    % Create a wavefunction with the desired energy
    phi = exp(-ka*(k-k0).^2).*exp(-1i*const*k.^2); 
    psi = ifft(phi);
    psi = psi/sqrt(psi'*psi*dx);
    % Change direction randomly
    if randi([0,1]) == 1
        psi = conj(psi);
    end
    
    % Create the time evolution operator
    % Fin-dif Laplacian and Hamiltonian
    e = ones(N,1); 
    Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
    H = -(1/2)*Lap + spdiags(U,0,N,N);
    
    % Time displacement operator E = exp(-iHdT/hbar)
    E = expm(-1i*full(H)*dt/hbar);
    
    if available
        electrons.psi(:, newID) = psi;
        electrons.evolve(:,:,newID) = E;
    else
        electrons.psi = [electrons.psi, psi];
        electrons.evolve = cat(3, electrons.evolve, E);
    end
    
    %avgE = norm(phi'*0.5*diag(k.^2,0)*phi*dk/(phi'*phi*dk))
end

%% Get new photon
% Uses AM1.5G spectrum as a probability density to pull a new photon energy
function [energy] = getRandomEnergyPhoton(wavelength, irradiance, hbar)
    wavelength = randsample(wavelength, 1, true,  irradiance)*1e-9;
    energy = hbar*(299792458)/wavelength;
end

%% Generate Potential
% Generates the potential energy function representing the given data
function [U,structureSize] = generatePotential(energy, voltage, N, q)
    % Set up the matrix
    U = zeros(1, N);
    
    % Get the size in nm of the structure;
    num = size(energy);
    num = num(1); % number of films
    structureSize = 0;
    for i = 1:num
        structureSize = structureSize + energy(i,2);
    end
    
    % Divide up the potential into the regions
    currentDivider = 1;
    for i = 1:num
        nextDivider = floor(N*(energy(i,2) / structureSize))+currentDivider;
        if i == num
            nextDivider = N;
        end
        U(1,currentDivider:nextDivider) = energy(i,1);
        currentDivider = nextDivider;
    end
    
    % Add external voltage
    for i = 1:N
       U(1,i) = U(1,i) + voltage * (i/N);
    end
    
    U = U.';
end

%% Utility plot function
% Plots functions with independant y-axes but shared x-axes in the same
% MATLAB figure. 
function [] = multiPlot(figureNum,x1, y1, x2, y2, ymax)
    figure(figureNum);
    cla reset
    line(x1, y1, 'Color', 'r');
    haxes1 = gca;
    set(haxes1, 'xtick', [], ...
        'xticklabel', [], ...
        'ytick', [], ...
        'yticklabel', [], ...
        'YColor', 'r');
    ylim([0, ymax]);
    
    
    haxes1_pos = get(haxes1,'Position'); % store position of first axes
    haxes2 = axes('Position',haxes1_pos,...
              'XAxisLocation','top',...
              'YAxisLocation','right',...
              'Color','none');
          
    line(x2,y2,'Parent',haxes2,'Color','k')
end

%% Normalize
% Normalizes a given discretized function
function psi = normalize(psi, dx)
    summa = dx * sum(psi.*conj(psi)); 
    weight = sqrt(1/summa);
    psi = psi * weight;
end
