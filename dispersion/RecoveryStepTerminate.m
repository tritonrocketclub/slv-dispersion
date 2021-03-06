function [ cont ] = RecoveryStepTerminate( sv, varargin )
% RECOVERYSTEPTERMINATE 
% 
% Objective: Terminate the recovery step function
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

cont = sv(3) > 0;

end

