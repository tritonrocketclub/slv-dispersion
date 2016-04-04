function [ newsv ] = RK45Chute( sv, args )
% RK45Step
% 
% Objective: Calculate an individual step of 4th order RK method to
%   solve an ODE representing the 3DOF equations of motion of
%   a rocket in the recovery phase
%
% input variables:
%   sv - row vector, name stands for State Vector
%       contains position, velocity, and time in the format
%       [x,y,z,vx,vy,vz,t]
%   args - cell array, variable input arguments to pass to the function
%       {1} - must be a function handle that calculates acceleration
%       Note that this is NOT args as args was already passed to a
%       previous function, and this would cause nested cell arrays
%
% output variables:
%   newsv - row vector, the new state vector for the fragment
%
% functions called:
%   none
%

%
% Initialize global structures
%
global rocketProp;
global simuProp;
%
% Consolidate state vector terms
%
pos = [sv(1), sv(2), sv(3)];
vel = [sv(4), sv(5), sv(6)];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate RK45
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate an intelligent time step based on the magnitude
% of the acceleration
% Bound the step between 1/1000 and 1/2
% Create h/2 variable
%
h = simuProp.recoveryRes;
if h < 10^-3
    h = 10^-3;
elseif h > .5
    h = .5;
end
h2 = h/2;
%
% Calculate position K1 term
%
vK1 = ChuteVel( sv, args );
pK1 = h.*vK1;
%
% Calculate position K2 term
% 
svK2 = sv + [pK1./2, vK1./2, h2];
vK2 = ChuteVel( svK2, args );
pK2 = h.*vK2;
%
% Calculate position K3 term
% 
svK3 = sv + [pK2./2, vK2./2, h2];
vK3 = ChuteVel( svK3, args );
pK3 = h.*vK3;%
%
% Calculate position K4 term
% 
svK4 = sv + [pK3, vK3, h];
vK4 = ChuteVel( svK4, args );
pK4 = h.*vK4;
%
% Calculate final step position
%
vstep = (1/6).*(vK1 + 2.*(vK2 + vK3) + vK4);
pstep = pos + (1/6).*(pK1 + 2.*(pK2 + pK3) + pK4);
%
% Format output state vector
%
t = sv(7) + h;
newsv = [pstep, vstep, t];