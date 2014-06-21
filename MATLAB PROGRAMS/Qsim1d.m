function [] = Qsim1d()
% Yee Algorithm (ie, the first thing anyone thinks of trying)
% 1. Replace all the derivatives in Ampere’s and Faraday’s laws with finite 
% differences. Discretize space and time so that the electric and magnetic 
% fields are staggered in both space and time.
%
% 2. Solve the resulting difference equations to obtain “update equations” 
% that express the (unknown) future fields in terms of (known) past fields.
%
% 3. Evaluate the magnetic fields one time-step into the future so they are 
% now known (effectively they become past fields).
%
% 4. Evaluate the electric fields one time-step into the future so they are
% now known (effectively they become past fields).
%
% 5. Repeat the previous two steps until the fields have been obtained over 
% the desired duration.


vidname = strcat('C:\Users\Brandon\Desktop\',strrep(datestr(clock()),':',''));
aviobj = avifile(vidname, 'fps', 60);

L = 1;
N = 500;
iterations = 1000;
frame_rate = 60;
x = linspace(0, L, N);
m = 9.11e-31;
dx = x(2) - x(1);

% Gaussian initial wf
psi = gaussmf(x, [L/100, L/2]).* exp(-(1i*100*x));

norm = sum(psi.*conj(psi));
psi = psi / sqrt(norm);
% Potential
v = x*0;

% Constants
epsilon_relative = 1;
mu_relative = 1;
hbar = 1.054571726e-34;
mu0 = 4 * pi * 1e-7; 
e0 = 8.854e-12;
n0 = sqrt(mu0/e0); % Characteristic impedance of free space
c = 1/sqrt(mu0*e0);

% Just say you cant move a space faster than a little less than c, I guess
dt = 1000 * dx^2;

% "spatial admittance"
imps = 1i * hbar * dt / 2 / m / dx^2;
% "potential admittance"
impv = 1i * dt / hbar;

fig_h = figure('Visible', 'off');

% Fin-diff
for t = 1:iterations
    
    psi_diff = x*0;
    for m = 1:N
        if m ~= 1 && m ~= N
            psi_diff(m) = imps *(psi(m+1)-2*psi(m)+psi(m-1))-impv*v(m)*psi(m);
        elseif m == 1
            psi_diff(m) = imps *(psi(m+1)-2*psi(m)+0)-impv*v(m)*psi(m);
        elseif m == N
            psi_diff(m) = imps *(0+2*psi(m)+psi(m-1))-impv*v(m)*psi(m);
        end
    end
    
    psi = psi + psi_diff;
    plot(x, real(psi));
%     xlim([0,L]);
%     ylim([-5,5]);
    aviobj = addframe(aviobj, fig_h);
    
end
aviobj = close(aviobj)

end

