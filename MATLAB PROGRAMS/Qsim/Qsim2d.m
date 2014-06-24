function [] = Qsim2d()
% A 2D version of Qsim1d

L = 100;
N = 1000;
iterations = 5000;
x = linspace(0, L, N)';
% m = 9.11e-31;
m = 1;
dx = x(2) - x(1);

% Gaussian initial wf 
omega = 15*pi;
psi = (gaussmf(x, [L/100, L/4])).*exp(-1i*omega.*x);
psi = psi / sqrt(psi'*psi*dx);

% Potential
%U = -(1.602e-19)*gaussmf(x, [L/50, L/2]);
%U = x*0;
%U = (x-L/2).^2;
U = 100*heaviside(x-L/2);
% %U(floor(N/2):end) = (1.602e-19);

% Constants
epsilon_relative = 1;
mu_relative = 1;
% hbar = 1.054571726e-34;
hbar = 1;
mu0 = 4 * pi * 1e-7; 
e0 = 8.854e-12;
n0 = sqrt(mu0/e0); % Characteristic impedance of free space
c = 1/sqrt(mu0*e0);

% Just say you cant move a space faster than a little less than c, I guess
dt =  0.01;

% This actually makes sense now
e = ones(N,1);
Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
H = -((hbar^2) / 2 / m) * Lap + spdiags(U,0,N,N);

% Exponential matrix formulation of the time evolution operator
E = expm(-1i * full(H) * dt / hbar);

% Fin-diff
for t = 1:iterations
    
    psi = E*psi;
    
    clc;
    fprintf('Percent complete: %5.2f \n',100 * t / iterations);
    fprintf('Transmission percent: %5.2f \n', 100 * sum(psi(floor(N/2):end).*conj(psi(floor(N/2):end)))/sum(psi.*conj(psi)));
    fprintf('Reflection percent: %5.2f \n', 100 * sum(psi(1:floor(N/2)).*conj(psi(1:floor(N/2))))/sum(psi.*conj(psi)));
    figure(1)
    plot(x, psi.*conj(psi));
    
end


end

