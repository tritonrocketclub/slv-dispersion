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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate K1
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculate acceleration
%
aK1 = BreakupAccel(sv, fmass, Cd, A, simp, buf);
%
% Calculate an intelligent time step based on the magnitude
% of the acceleration
% Bound the step between 1/1000 and 1/2
% Create h/2 variable
%
res = simp.res;
h = 100/(res*sum(aK1(1)^2 + aK1(2)^2 + aK1(3)^2));
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
svK2 = sv + [vK1./2, pK1./2, h2];
%
% Calculate acceleration with K2 state vector
%
aK2 = BreakupAccel(svK2, fmass, Cd, A, simp, buf);
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
svK3 = sv + [vK2./2, pK2./2, h2];
%
% Calculate acceleration with K3 state vector
%
aK3 = BreakupAccel(svK3, fmass, Cd, A, simp, buf);
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
svK4 = sv + [vK3, pK3, h];
%
% Calculate acceleration with K4 state vector
%
aK4 = BreakupAccel(svK4, fmass, Cd, A, simp, buf);
%
% Calculate velocity K4 term
%
vK4 = h.*aK4;
vel4 = [svK4(4), svK4(5), svK4(6)];
%
% Calculate position K4 term
%
pK4 = h.*vel4;
%
% Calculate the sum RK45 terms
%
vstep = vel1 + (1/6).*(vK1 + 2.*(vK2 + vK3) + vK4);
pstep = pos1 + (1/6).*(pK1 + 2.*(pK2 + pK3) + pK4);
%
% Format output state vector
%
t = sv(7) + h;
newsv = [pstep, vstep, t];
end