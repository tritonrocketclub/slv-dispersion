function [ newsv ] = RK45Wind( sv, args )
% RK45Step
% 
% Objective: Calculate an individual step of 4th order RK method to
%   solve an ODE representing the 3DOF equations of motion of
%   a rocket with offset due to wind action
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
% Get the acceleration function
%
fn = args{1};
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K1
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculate acceleration
%
aK1 = fn(sv, args);
%
% Calculate an intelligent time step based on the magnitude
% of the acceleration
% Bound the step between 1/1000 and 1/2
% Create h/2 variable
%
h = 100/(simuProp.resolution*sum(aK1(1)^2 + aK1(2)^2 + aK1(3)^2));
if h < 10^-3
    h = 10^-3;
elseif h > .5
    h = .5;
end
h2 = h/2;

%
% Calculate velocity K1 term
%
vK1 = h.*aK1;
%
% Calculate position K1 term
%
vel1 = [sv(4), sv(5), sv(6)];
pos1 = [sv(1), sv(2), sv(3)];
pK1 = h.*vel1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K2
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create new state vector for K2
%
svK2 = sv + [pK1./2, vK1./2, h2];
%
% Calculate acceleration with K2 state vector
%
aK2 = fn(svK2, args);
%
% Calculate velocity K2 term
%
vK2 = h.*aK2;
%
% Calculate position K2 term
%
vel2 = [svK2(4), svK2(5), svK2(6)];
pK2 = h.*vel2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K3
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create new state vector for K3
%
svK3 = sv + [pK2./2, vK2./2, h2];
%
% Calculate acceleration with K3 state vector
%
aK3 = fn(svK3, args);
%
% Calculate velocity K3 term
%
vK3 = h.*aK3;
vel3 = [svK3(4), svK3(5), svK3(6)];
%
% Calculate position K3 term
%
pK3 = h.*vel3;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K4
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create new state vector for K4
%
svK4 = sv + [pK3, vK3, h];
%
% Calculate acceleration with K4 state vector
%
aK4 = fn(svK4, args);
%
% Calculate velocity K4 term
%
vK4 = h.*aK4;
vel4 = [svK4(4), svK4(5), svK4(6)];
%
% Calculate position K4 term
%
pK4 = h.*vel4;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Sum Terms
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the wind at altitude
%
wvelx = interp1(simuProp.windsr, simuProp.windse, sv(3));
wvely = interp1(simuProp.windsr, simuProp.windsn, sv(3));
wvelz = 0;
wvel = [wvelx, wvely, wvelz];
%
% Calculate the sum RK45 terms
%
vstep = vel1 + (1/6).*(vK1 + 2.*(vK2 + vK3) + vK4);
pstep = pos1 + (1/6).*(pK1 + 2.*(pK2 + pK3) + pK4) + h.*wvel;
%
% Format output state vector
%
t = sv(7) + h;
newsv = [pstep, vstep, t];
end

