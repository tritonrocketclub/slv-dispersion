function [ output_args ] = LogProp( propertyName, propertyVal )
% LogProp
% 
% Objective: Log the latest property - index should match that of timestamp
%
% input variables:
%   propertyName - name of the property to log
%   propertyVal - value of the property 
%
% output variables:
%   sm - matrix, a matrix describing the rocket's trajectory with
%   multiple rows of state vectors
%
% functions called:
%   none

global simuProp
global logProp

i = logProp.([propertyName, '_int']);
logProp.(propertyName)(i) = propertyVal;
i = i + 1;
logProp.([propertyName, '_int']) = i;

end