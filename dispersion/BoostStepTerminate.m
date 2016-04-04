function [ cont ] = BoostStepTerminate( sv, varargin )
% BOOSTSTEPTERMINATE 
% 
% Objective: Terminate the boost step function
%
% input variables:
%   sv - row vector, name stands for State Vector
%       contains position, velocity, and time in the format
%       [x,y,z,vx,vy,vz,t]
%   varargin - cell array, variable input arguments to pass to the function
%       {1} - will be a function handle that calculates acceleration,
%       neglected
%       {2} - will be the flight path angle
%
% output variables:
%   cont- whether or not to continue calculating steps
%
% functions called:
%   none
%

%
% Initialize global structures
%
global rocketProp
cont = sv(7) < rocketProp.burnTime;

end

