function [ outa ] = BreakupAccel( sv, fmass, Cd, A, simp, buf )
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
% Apply RungeKutta45
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K1
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the drag force
%
pos = [sv(1), sv(2), sv(3)];
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
% Sum the accelerations
%
ax = fD(1)/fmass + fBu(1)/fmass;
ay = fD(2)/fmass + fBu(2)/fmass;
az = fD(3)/fmass + fBu(3)/fmass - ig;
%
% Create output acceleration vector
%
outa = [ax, ay, az];