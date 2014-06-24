function [] = Qsim2d()
% A 2D version of Qsim1d
% Simply

% The sim will be 
L = 100;
N = 250;
iterations = 5000;

x = linspace(0,L,N)';
y = linspace(0,L,N)';
% m = 9.11e-31;
m = 1;
a = x(2) - x(1);

% Gaussian initial wf 
omega = 16*pi;
ymod = gaussmf(y, [L/100, L/2]);
for n = 1:N
    psi(:,n) = gaussmf(x, [L/100, L/4]) * ymod(n) .*exp(+1i*omega.*x);
end

% Normalize
psi = psi / sqrt(norm2(psi,a));

% Potential
U = zeros(N,N);
U(floor(N/2):end,:) = 1;
figure(2); surf(U);


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


% Fin-diff
for t = 1:iterations
    
    psi = psi + (dt * 1i * hbar / 2 /m) * laplacian(psi, a) - (dt * 1i / hbar) * U * psi ;
    psi = psi / sqrt(norm2(psi,a));
    clc;
    fprintf('Percent complete: %5.2f \n',100 * t / iterations);
%     fprintf('Transmission percent: %5.2f \n', 100 * sum(psi(floor(N/2):end).*conj(psi(floor(N/2):end)))/sum(psi.*conj(psi)));
%     fprintf('Reflection percent: %5.2f \n', 100 * sum(psi(1:floor(N/2)).*conj(psi(1:floor(N/2))))/sum(psi.*conj(psi)));
    figure(1)
    h = surf(psi.*conj(psi));
    set(h, 'edgecolor', 'none');
end


end

function out = norm2(in, a)

out = sum(sum((in.*conj(in)*a^2)));

end

function Lap = laplacian(in,dx)
% Computes the 2d laplacian.
% Annoyingly, matricies in matlab don't seem to generalize up to nD so we
% cant use that wonderful matrix exponential method

[M,N] = size(in);
Lap = zeros(M,N);

for m = 1:M
    for n = 1:M
        if (m ~= M)
            mp1 = in(m+1,n);
        else
            mp1 = 0;
        end
        if (m ~= 1)
            mn1 = in(m-1,n);
        else
            mn1 = 0;
        end
        if (n ~= N)
            np1 = in(m,n+1);
        else
            np1 = 0;
        end
        if (n ~= 1)
            nn1 =  in(m,n-1);
        else
            nn1 = 0;
        end
        Lap(m,n) = mp1 + mn1 + np1 + nn1 - 4*in(m,n);
    end
end
Lap = Lap / dx^2;

end

