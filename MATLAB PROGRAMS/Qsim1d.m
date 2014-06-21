function [] = Qsim1d()

L = 100; % Interval Length
N = 1500; % No of points
x = linspace(0,L,N)'; % Coordinate vector
dx = x(2) - x(1); % Coordinate step

% Parameters for making intial momentum space wavefunction phi(k)
ko = 2; % Peak momentum
a = 20; % Momentum width parameter
dk = 2*pi/L; % Momentum step
km=N*dk; % Momentum limit
k=linspace(0,+km,N)'; % Momentum vector

% Make psi(x,0) from Gaussian kspace wavefunction phi(k) using
% fast fourier transform :
phi = exp(-a*(k-ko).^2).*exp(-1i*10*k.^2); % unnormalized phi(k)
psi = ifft(phi); % multiplies phi by expikx and integrates vs. x
psi = psi/sqrt(psi'*psi*dx); % normalize the psi(x,0)
plot(x, psi); 
V = zeros(1,N);


iterations = 2000;
hbar = 1;
m = 1;
dt = dx^2;

for t = 1:iterations
    delta_psi = zeros(N,2);
    psi(1,1) = 0;
    psi(end, 1) = 0;
    % Loop through x and compute the change
    for n = 1:N
        if n ~= N && n~= 1
            delta_psi(n,1) = (-1i * dt * hbar / m) * ( psi(n+1) - 2*psi(n) +  psi(n-1))/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
        if n == N
            delta_psi(n,1) = (-1i * dt * hbar / m) * ( 0 - 2*psi(n) +  psi(n-1))/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
        if n == 1
            delta_psi(n,1) = (-1i * dt * hbar / m) * ( psi(n+1) - 2*psi(n) +  0)/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
    end
    for n = N:1
        if n ~= N && n~= 1
            delta_psi(n,2) = (-1i * dt * hbar / m) * ( psi(n+1) - 2*psi(n) +  psi(n-1))/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
        if n == N
            delta_psi(n,2) = (-1i * dt * hbar / m) * ( 0 - 2*psi(n) +  psi(n-1))/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
        if n == 1
            delta_psi(n,2) = (-1i * dt * hbar / m) * ( psi(n+1) - 2*psi(n) +  0)/dx + (-1i * dt / hbar) * V(n) * psi(n); 
        end
    end
    delta_psi = (delta_psi(:,1) + delta_psi(:,2)) / 2;
    psi = psi + delta_psi;
    norm = sum(psi.*conj(psi));
    %psi =  psi / norm;
    figure(1);
    cla reset
    plot(x, psi.*conj(psi));
end

end

