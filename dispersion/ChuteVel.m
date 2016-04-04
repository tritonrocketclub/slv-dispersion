function [ vv ] = ChuteVel( sv, args )
% COASTACCEL
% 
% Objective: Calculate the velocity for an individual term in the RK
%   calculation of recovery/chute velocity
%
% input variables:
%   sv - row vector, name stands for State Vector
%       contains position, velocity, and time in the format
%       [x,y,z,vx,vy,vz,t]
%   args - cell array, variable input arguments to pass to the function
%       {1} - will be a function handle that calculates acceleration,
%       neglected
%       {2} - will be the flight path angle
%
% output variables:
%   vv- row vector, contains 3 velocity terms [vx, vy, vz]
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
% Organize rocket terms
%
pos = [sv(1), sv(2), sv(3)];
vel = [sv(4), sv(5), sv(6)];
imass = rocketProp.dryMass;
%
% Get the current gravity, air density, and temperature
%
ig = interp1(simuProp.gravityr, simuProp.gravity, sv(3));
irho = interp1(simuProp.rhor, simuProp.rho, sv(3));
itemp = interp1(simuProp.tempr, simuProp.temp, sv(3));
itemp = itemp + 273.15; %convert *C to K
%
% Calculate the wind velocity
%
wvelx = interp1(simuProp.windsr, simuProp.windse, sv(3));
wvely = interp1(simuProp.windsr, simuProp.windsn, sv(3));
%
% Calculate the terminal velocity with the parachute
%
Cd = rocketProp.chuteDrag;
R = rocketProp.chuteRadius;
r = rocketProp.chuteSpillRadius;
cvel = sqrt(2*imass*ig/(irho*pi*(Cd*R^2*(cos((r*pi)/(2*R)))^2)));
%
% Export velocity vector
%
vv = [wvelx, wvely, -cvel];

end