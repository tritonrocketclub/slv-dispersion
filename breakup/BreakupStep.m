function [ newsv ] = BreakupStep( sv, fmass, Cd, A, simp, buf )
% BREAKUPSTEP
% 
% Objective: Calculate an individual step of 2nd order RK method to
%   solve an ODE representing the 3DOF equations of motion of
%   a fragment of a rocket fragment
%
% input variables:
%   sv - row vector, name stands for State Vector
%       contains position, velocity, and time in the format
%       [x,y,z,vx,vy,vz,t]
%   fmass - number, fragment mass
%       Used to calculate equations of motion
%   Cd - number, the drag coefficient for the fragment
%   A - number, the reference area for the fragment
%   simp - struct, contains simulation properties
%   buf - struct, contains break up force properties
%
% output variables:
%   newsv - row vector, the new state vector for the fragment
%
% functions called:
%   none
%
res = simp.res;
gravity = simp.gravity;
rho = simp.rho;
%
% Calculate current gravity and air density
%
ig = interp1(gravity(:,1), gravity(:,2), sv(3)); %[m/s^2] Precise acceleration due to gravity
irho = interp1(rho(:,1), rho(:,2), sv(3));
%
% Calculate the drag force
%
vel = [sv(4), sv(5), sv(6)];
fmD = Cd*irho*A/2*(vel(1)^2+vel(2)^2+vel(3)^2);
%
% Get unit vector in direction of velocity
%
nvel = norm(vel);
if nvel == 0
    vDir = [0, 0, 0];
else
    vDir = vel/norm(vel);
end
%
% Calculate drag force in opposite direction to motion
%
fD = fmD*-vDir;
%
% Calculate the gravitational force
%
fmG = ig*fmass;
fG = [0,0,-fmG];
%
% Calculate the explosive breakup force using inverse square
%
busv = buf.busv;
fmBu = buf.force;
dirBu = buf.direction;
distsqr = 1+(sv(1)-busv(1))^2 + (sv(2)-busv(2))^2 + (sv(3)-busv(3))^2;
fBu = dirBu*fmBu./distsqr;
%
% Apply 2nd Order Runge Kutta Method to solve ODE
%
% Acceleration
ax = fD(1)/fmass + fBu(1)/fmass;
ay = fD(2)/fmass + fBu(2)/fmass;
az = fD(3)/fmass + fBu(3)/fmass - ig;
%
% Calculate an intelligent time step based on the magnitude
% of the acceleration
% Bound the step between 1/1000 and 1/2
%
h = 100/(res*sum(ax^2 + ay^2 + az^2));
if h < 10^-3
    h = 10^-3;
elseif h > .5
    h = .5;
end
%
% Calculate h/2
%
h2 = h/2;
% Velocity
vx = vel(1) + h/2*ax;
vy = vel(2) + h/2*ay;
vz = vel(3) + h/2*az;
% Position
px = sv(1) + h/2*vx;
py = sv(2) + h/2*vy;
pz = sv(3) + h/2*vz;
%
% Calculate drag force again using new velocity
% Make sure to recalculate direction of force
% [solution to an obscure bug]
%
fmD = Cd*irho*A/2*(vx^2+vy^2+vz^2);
vel = [vx, vy, vz];
nvel = norm(vel);
if nvel == 0
    vDir = [0, 0, 0];
else
    vDir = vel/norm(vel);
end
fD = fmD*-vDir;
%
% Calculate breakup force again using new position
%
distsqr = 1 + (px-busv(1))^2 + (py-busv(2))^2 + (pz-busv(3))^2;
fBu = dirBu*fmBu./distsqr;
% decay = exp(-(sv(7)+h/2)/busv(7))*(1-(sv(7)+h/2)/busv(7));
% fBu = fBu*decay;
%
% Recalculate ODE
%
% Acceleration
ax = fD(1)/fmass + fBu(1)/fmass;
ay = fD(2)/fmass + fBu(2)/fmass;
az = fD(3)/fmass + fBu(3)/fmass - ig;
% Velocity
vx = vel(1) + h*ax;
vy = vel(2) + h*ay;
vz = vel(3) + h*az;
% Position
px = sv(1) + h*vx;
py = sv(2) + h*vy;
pz = sv(3) + h*vz;
% Time
t = sv(7) + h;
%
% Format the output state vector
%
newsv = [px, py, pz, vx, vy, vz, t];
%
% End of function BreakupStep.m
%
end
