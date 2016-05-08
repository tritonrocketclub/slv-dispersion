function [ sm ] = RocketTrajectory( sv, useFPA )
% BOOST
% 
% Objective: Take an initial condition state vector and calculate the
%   total trajectory of the rocket
%
% input variables:
%   sv - row vector, a state vector describing the rocket's initial
%   position
%   useFPA - bool, whether or not to use the FPA or default to zero
%
% output variables:
%   sm - matrix, a matrix describing the rocket's trajectory with
%   multiple rows of state vectors
%
% functions called:
%   none

%
% Establish global variables
%
global rocketProp;
global simuProp;
%
% Calculate flight path angle vector vertical
%
if useFPA
    gamma = GetSTDFPA();
    FPAv = normrnd(0, gamma);
    max(min(FPAv, pi/2),-pi/2);
    %
    % General FPA lateral
    %
    FPAl = 2*pi*rand();
    %
    % Convert to cartesian coords
    %
    x = cos(FPAl)*sin(FPAv);
    y = sin(FPAl)*sin(FPAv);
    z = cos(FPAv);
    FPA = [x,y,z];
%
% Reset to vertical if not using FPA
%
else
    FPA = [0,0,1];
end
%
% Simulate the boost phase
%
fprintf('Calculating boost phase...');
tstamp = cputime;
boostm = IntegrateStepFunction( @RK45Wind, @BoostStepTerminate, sv, @BoostAccel, FPA );
simuProp.burnOut = boostm(end,:); 
fprintf('done! Took %0.3fs\n', cputime-tstamp);
%
% Simulate the coast phase
%
fprintf('Calculating coast phase...');
tstamp = cputime;
sv = boostm(end,:);
coastm = IntegrateStepFunction( @RK45Wind, @CoastStepTerminate, sv, @CoastAccel);
fprintf('done! Took %0.3fs\n', cputime-tstamp);
%
% Simulate the recovery phase
%
fprintf('Calculating recovery phase...');
tstamp = cputime;
sv = coastm(end,:);
recoverym = IntegrateStepFunction( @RK45Chute, @RecoveryStepTerminate, sv);
fprintf('done! Took %0.3fs\n', cputime-tstamp);
%
% Format output matrix
%
sm = [boostm;coastm;recoverym];
end