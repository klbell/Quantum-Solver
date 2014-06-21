
clc; clear; % clear console and variables

SIZE = 200; % array size

ez = zeros(SIZE,1); % Empty array for electric field
hy = ez; % Empty array for magnetic field
imp0 = 377.0; % Characteristic impedance of free space

maxTime = 500; % Time steps to evaluate

for qTime = 1:maxTime % from 1 to maximum time steps
    
   for mm = 1:SIZE-1 % for the size of array     
       hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) / imp0; % Update mag field      
   end
   
   for mm = 2:SIZE % for the size of array
       ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * imp0; % Update e field
   end 
   
    ez(50) = ez(50) +( exp(-(qTime - 30) * (qTime - 30) / 100)); % initial e field 
    
    
    plot(1:200,ez) % plot result   
    axis([0 200 -2 2]);
    pause(1/60); % space plots out for display
    
end