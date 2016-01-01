function [ trajectories ] = DispersionAnalysis( varargin )
% DISPERSIONANALYSIS
% 
% Objective: Calculate an individual step of Euler's method to
%   solve an ODE representing the 3DOF equations of motion of
%   a fragment of a rocket fragment
%
% variable input arguments:
%   {1} - number, resolution; default is 100
%   {2} - number, number of flights to simulate
%   {3} - number, number of fragments to simulate per flight
%
% output variables:
%   trajectories - cell array, contains all state vectors for all
%       breakup trajectories
%
% functions called:
%   BreakupTrajectory
%   BreakupStep
%
% Set the random number generator to the default seed
% This is for repeatability. It will not negatively impact
% the randomness of the function or the usefulness of the data
%
rng default;
%
% Handle input variables
%
if nargin == 3
    res = varargin{1};
    numsims = varargin{2};
    numfrag = varargin{3};
else
    res = input('Specify resolution (default 5): ');
    numsims = input('Specify number of iterations: ');
    numfrag = input('Specify number of fragments: ');
end
%
%  Load atmospheric data
%
fprintf('Loading atmospheric data...');
fileID = fopen('StandardAtmosphere.txt');
errmsg = 'No StandardAtmosphere.txt found in local folder. Please specify a new file.';
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end     
end
% Simple file format columns with %-style comments
AtmProp = textscan(fileID, '%f\t%f\t%f\t%f\t%f\t%f', 'CommentStyle', '%');
fclose(fileID);
fprintf('done!\n');
%
% Load rocket flight path data
%
fprintf('Loading flight path data...');
fileID = fopen('FlightData.txt');
errmsg = 'No FlightSim.csv found in local folder. Please specify a new file.';
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end     
end
% Simple file format columns with #-style comments
FlightSim = textscan(fileID, '%f\t%f\t%f\t%f\t%f\t%f', 'CommentStyle', '#');
fclose(fileID);
fprintf('done!\n');
%
% Format simulation properties structure
% Use atmospheric data and resolution input
%
gravity = [AtmProp{1}, AtmProp{3}];
rho = [AtmProp{1}, AtmProp{5}];
pressure = [AtmProp{1}, AtmProp{4}];
simp = struct('res', res, 'gravity', gravity, 'rho', rho, 'pressure', pressure);
%
% Format flight simulation data
%
time = FlightSim{1};
alt = FlightSim{2};
ldist = FlightSim{4};
ldir = FlightSim{5};
vvel = FlightSim{3};
lvel = FlightSim{6};
% Convert lateral distance and velocity to cartesian coords
xdist = ldist.*cos(ldir);
ydist = ldist.*sin(ldir);
xvel = lvel.*cos(ldir);
yvel = lvel.*sin(ldir);
% Convert flight data to state vector format
flightdata = [xdist, ydist, alt, xvel, yvel, vvel, time];
%
% Plot the successful rocket flight path
%
plot3( xdist, ydist, alt, 'LineWidth', 10, 'color', 'g');
hold on;
%
% Pre-allocate trajectory cell array
%
trajectories = cell(1, numsims*numfrag);
%
% Create fragment data text file
%
fragfileID = fopen('FragmentData.txt', 'w');
%
% Store cpu time for calculating total run time
% Print to console
%
tasktime = cputime;
fprintf('Calculating trajectories...');
%
% Calculate the minimum and maximum indices for
% the flight data, and the maximum noise size
%
minflight = 10;
maxflight = size(flightdata, 1);
flightnoise = (maxflight-minflight)/numfrag;
%
% Do the trajectory calculations for j iterations
% with i fragments and k total trajectories.
%
trajectories = {};
k = numsims*numfrag;
for j = 1:numsims
    flightpts = linspace(minflight, maxflight, numfrag);
    jtrajectories = {};
    parfor i = 1:numfrag
        % Add random noise to each flight point
        % This way the BUSVs vary for each iteration
        flightpts(i) = flightpts(i) + flightnoise*(2*rand()-1);
        l = round(flightpts(i));
        if l > maxflight
            l = maxflight;
        elseif l < minflight
            l = minflight;
        end
        busv = flightdata(l, :);
        states = BreakupTrajectory( busv, 40, simp, fragfileID, i+(j-1)*numfrag );
        statex = abs(states(end, 1));
        statey = abs(states(end, 2));
        statez = abs(states(end, 3));
        if statex < 10^4 && statey < 10^4 && statez < 10^5
            plot3(states(:,1), states(:,2), states(:,3));
        else
            for m = size(states, 1):-1:1
                if (states(m, 1) < 10^4) && (states(m, 2) < 10^4) && (states(m, 2) < 10^6)
                    plot3(states(1:res:m,1), states(1:res:m,2), states(1:res:m,3));
                    break
                end
            end
        end
        jtrajectories{i}=states;
        hold on;
    end
    trajectories = [trajectories,jtrajectories];
    fprintf('\nCompleted %d iterations at %2.2f seconds...', j, cputime-tasktime);
end
%
% Close the fragment data file
%
fclose(fragfileID);
%
% Clean up the plot axes
%
lim = axis;
axis( [lim(1), lim(2), lim(3), lim(4), 0, lim(6)] );
%
% Calculate elapsed time and print to console
%
fprintf('done!\n');
tasktime = cputime-tasktime;
strtime = sprintf('Total elapsed time: %2.2fs\n', tasktime);
fprintf(strtime);
%
% Calculate mean distance and standard deviation
% Print these to console
%
dist = zeros(1, k);
allx = zeros(1, k);
ally = zeros(1, k);
for i = 1:k
    allx(i) = trajectories{i}(end, 1);
    ally(i) = trajectories{i}(end, 2);
end
meanx = mean(allx);
meany = mean(ally);
stdx = std(allx);
stdy = std(ally);
for i = 1:k
    dist(i) = sqrt((allx(i) - meanx)^2 + (ally(i) - meany)^2);
end
distmean = mean(dist);
diststd = std(dist);
fprintf('Mean x: %2.2fm\n', meanx);
fprintf('Mean y: %2.2fm\n', meany);
fprintf('Mean distance: %2.2fm\n', distmean);
fprintf('StdDev on distance: %2.2fm\n', diststd);
fprintf('3sigma distance: %2.2fm\n', 3*diststd+distmean);
%
% Plot circles
%
t = linspace(0, 2*pi, 200);
z = zeros(1, 200);
plot3( meanx + diststd.*cos(t), meany + diststd.*sin(t), z, 'Color', 'c', 'LineStyle', '--');
plot3( meanx + 2.*diststd.*cos(t), meany + 2.*diststd.*sin(t), z, 'Color', 'c', 'LineStyle', '--');
plot3( meanx + 3.*diststd.*cos(t), meany + 3.*diststd.*sin(t), z, 'Color', 'c', 'LineStyle', '--');
%
% Inform user that dispersion analysis is completed
%
str1 = 'Dispersion simulation completed.';
h = msgbox({'Dispersion simulation completed.' strtime});
%
% End of function DispersionAnalysis.m
%
end