function [ texdata ] = ReadFile( infile, easyname, cols )
% READFILE
% 
% Objective: Open a file and output a matrix containing the data
%   for processing
%
% input variables:
%   infile - string, the system filename from which to read, a txt file
%   easyname - string, the name to display when opening the file
%   cols - integer, number of columns to read. 0 if "name = var" format
%
% output variables:
%   texdata - a matrix containing the read data
%
% functions called:
%   none
%

%
% Print messages
%
fprintf('Loading %s...', easyname);
fileID = fopen(infile);
errmsg = sprintf('No %s found in local folder. Please specify a new file.', infile);
%
% Check that the default file opened properly
% While no valid file is open, ask for a new one
% If the file we asked for is invalid, try appending .txt
%
while fileID < 0
    disp(errmsg);
    filename = input('Open file: ', 's');
    [fileID, errmsg] = fopen(filename);
    if fileID == -1
        filename = sprintf('%s.txt', filename);
        [fileID, errmsg] = fopen(filename);
    end
    fprintf('Loading %s...', easyname);
end
%
% Scan the properties text file for contents, string = value
% Use comments styled like Matlab, %
%
if (cols == 0)
    readstr = '%s = %f';
else
    readstr = '%f';
    if (cols > 1)
        for i = 2:cols
            readstr = [readstr, '\t%f'];
        end
    end
end
texdata = textscan(fileID, readstr, 'CommentStyle', '%');
fclose(fileID);
fprintf('done!\n');
%
% End function
%
end

