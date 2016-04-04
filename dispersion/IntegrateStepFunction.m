function [ sm ] = IntegrateStepFunction( fn, fnt, sv, varargin )
% INTEGRATESTEPFUNCTION
% 
% Objective: Run the integration process for a generic step function
%
% input variables:
%   fn - function handle, the function that calculates a given step
%   fnt - the terminal condition, accepts sv and all varargin
%   sv - row vector, a state vector describing the rocket's initial
%       position
%   varargin - any arguments to pass to the step function
%
% output variables:
%   sm - matrix, a matrix describing the rocket's trajectory with
%       multiple rows of state vectors
%
% functions called:
%   none
%

%
% Initialize global structures
%
global simuProp
%
% Integrate until terminal condition
%
sm = zeros(simuProp.resolution,7); %Pre-allocate space for trajectory states
sm(1, :) = fn( sv, varargin );
i = 1;
k = 2;
while fnt(sv, varargin)
    % Pre-allocate additional space
    if i > simuProp.resolution
        i = 1;
        sm = [sm; zeros(simuProp.resolution, 7)];
    end
    % Calculate new state vector
    sv = fn( sv, varargin );
    % Store new state row vector in output state matrix
    sm(k, :) = sv;
    i = i+1;
    k = k+1;
end
%
% Loop backwards through state matrix and remove any unnecessary
% pre-allocated rows
%
for i = size(sm, 1):-1:1
    if sm(i, 7) ~= 0
        sm = sm(1:i,:);
        break
    end
end
for i = size(sm, 1):-1:1
    if any(isnan(sm(i,:)));
        sm = sm(1:i-1,:);
        break
    end
end

end

