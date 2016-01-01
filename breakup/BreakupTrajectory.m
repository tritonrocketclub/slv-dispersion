function [ states ] = BreakupTrajectory( busv, mass, simp, fileID, n )
% BREAKUPTRAJECTORY
% 
% Objective: Simulate the flight path of a piece of debris
%       by applying Euler's method to an ODE representing the
%       3DOF equations of motion of the rocket fragment
%
% input variables:
%   busv - row vector, name stands for Break Up State Vector
%       contains position, velocity, and time in the format
%       [x,y,z,vx,vy,vz,t]
%   mass - number, rocket mass at break up state.
%       Used to calculate energy of imparted break up velocity
%   simp - struct, contains simulation properties
%   fileID - file handle, points to a text file where fragment
%       data will be stored 
%   n - number, fragment ID
%
% output variables:
%   states - matrix, where each row is a state vector
%       for the breakup trajectory
%
% functions called:
%   BreakupStep - calculates each step of the ODE
%
%
% Calculate a random mass for a fragment
% ranging from 0 to total rocket mass
% Uniformly distributed
%
fmass = 10^-5 + rand(1)*mass; 
%
% Calculate a reference area with uniform distribution
%
rocketA = 0.31; %[m^2] Profile area of the rocket
A = .31*rand(1);
%
% Calculate chemical energy at breakup
% http://www.philsrockets.org.uk/physics.pdf
%
k = 1.38*10^-23; %[J/K] Boltzmann's Constant
T = 1366;% [K] temperature of reaction
CE = (5/2)*(1.38*10^-23)*T;
%
% Assume all chemical energy is converted to
% kinetic energy in the fragment, and find
% the square of velocity
%
fvmag2 = 2*CE/fmass;
fvmag = sqrt(fvmag2);
fvmag = 300;
%
% Calculate breakup pressure
% 
pcato = 13.789*10^6;
pressure = simp.pressure;
patm = interp1(pressure(:,1), pressure(:,2), busv(3));
pdif = pcato-patm;
f = pdif*A;
% t = 0.001;
% fvmag = f*t/fmass;
% display(A);
% display(fvmag);
%
% Create an independent trivariate normal distribution
%
mu = [0, 0, 0];
sigma = [1,1,1]./3;
r = mvnrnd(mu, sigma);
%
% Create a perturbation velocity.
%
vp = fvmag*r*0;
buf = struct('force', f, 'direction', r, 'busv', busv);
%
% Add the perturbation velocity to the BUSV
%
sv = busv + [0,0,0,vp,0];
%
% Calculate a random drag coefficient centered on 0.75
% with normal distribution where 3sigma = 0.75
%
Cd = randn(1)*0.75/3 + 0.75;
if Cd < 0.005
    Cd = 0.005;
end
%
% Write fragment data to file
%
fprintf(fileID,'Generated fragment [%d] with parameters:\n', n);
fprintf(fileID,'Mass = %2.2f\n', fmass);
fprintf(fileID,'Percent Mass = %2.2f%%\n', fmass/mass*100);
fprintf(fileID,'Area = %2.2f\n', A);
fprintf(fileID,'Cd = %2.2f\n', Cd);
fprintf(fileID,'Perturbation = %2.2f, %2.2f, %2.2f\n', r);
fprintf(fileID,'Break up force = %2.2f\n', buf.force);
fprintf(fileID,'BUSV = %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f\n\n', busv);
%
% Get simulation resolution
%
res = simp.res;
%
% Loop until the altitude is below zero
%
states = zeros(res,7); %Pre-allocate space for trajectory states
states(1, :) = sv;
i = 1;
k = 2;
while sv(3) > 0
    % Pre-allocate additional space
    if i > res
        i = 1;
        states = [states; zeros(res, 7)];
    end
    % Calculate new state vector
    sv = BreakupStep( sv, fmass, Cd, A, simp, buf );
    % Store new state row vector in output state matrix
    states(k, :) = sv;
    i = i+1;
    k = k+1;
end
%
% Loop backwards through state matrix and remove any unnecessary
% pre-allocated rows
%
for i = size(states, 1):-1:1
    if states(i, 7) ~= 0
        states = states(1:i,:);
        break
    end
end
for i = size(states, 1):-1:1
    if any(isnan(states(i,:)));
        states = states(1:i-1,:);
        break
    end
end
%
% End of function BreakupTrajectory.m
%
end
