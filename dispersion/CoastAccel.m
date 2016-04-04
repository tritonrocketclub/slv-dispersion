function [ av ] = CoastAccel( sv, args )
% COASTACCEL
% 
% Objective: Calculate the acceleration for an individual term in the RK
%   calculation of coast velocity
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
%   av- row vector, contains 3 acceleration terms [ax, ay, az]
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
vmag = norm(vel);
vdir = vel/vmag;
%
% Get the current gravity, air density, and temperature
%
ig = interp1(simuProp.gravityr, simuProp.gravity, sv(3));
irho = interp1(simuProp.rhor, simuProp.rho, sv(3));
itemp = interp1(simuProp.tempr, simuProp.temp, sv(3));
itemp = itemp + 273.15; %convert *C to K
%
% Calculate the speed of sound and mach number
%
ic = 20.05*sqrt(itemp);
imach = vmag/ic;
%
% Calculate the drag coefficient. Make sure it's bounded.
% 
if imach <= rocketProp.drag(end, 1)
    iCd = interp1(rocketProp.dragr, rocketProp.drag, imach);
else
    iCd = rocketProp.drag(end);
end
%
% Calculate the drag force
%
fDm = 0.5*irho*vmag^2*rocketProp.area*iCd;
fDv = -1*vdir*fDm;
%
% Get the current mass from time
%
imass = rocketProp.dryMass;
%
% Calculate the gravitational force
%
fGm = ig*imass;
fGv = [0, 0, -fGm];
%
% Sum the accelerations
%
ax = fDv(1)/imass;
ay = fDv(2)/imass;
az = fDv(3)/imass - ig;
av = [ax, ay, az];
LogProp('acc', av(3));
LogProp('drag', fDv(3));
LogProp('dm', fDm);
LogProp('vel', vel(3));
LogProp('rho', irho)
LogProp('grav', fGm);
