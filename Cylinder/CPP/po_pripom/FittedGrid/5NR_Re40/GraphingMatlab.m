% Open the file
filename = 'flowdata_cylinder_NR25.000000_Re40.000000_Iter150000.txt';
fileID = fopen(filename, 'r');

% Read the first line (simulation time)
time = str2double(fgetl(fileID));

% Initialize an empty cell array to hold the data
data = {};

% Read the file line by line
line = fgetl(fileID);
while ischar(line)
    % Split the line by the ']' character and remove the empty string at the end
    line = strrep(line, '[', ''); % Remove the opening brackets
    line = strrep(line, ']', ''); % Remove the closing brackets
    entries = strsplit(line);
    
    % Process each entry in the line
    row_data = [];
    for i = 1:length(entries)
        % Convert the string to a numerical array
        row_data = [row_data; str2num(entries{i})]; % Append the row data to the matrix
    end
    
    % Append the row to the data array
    data = [data; row_data];
    
    % Read the next line
    line = fgetl(fileID);
end

% Close the file
fclose(fileID);

% Convert the cell array to a 2D matrix (each element is a row of 6 numbers)
data_matrix = cell2mat(data);

% Display the matrix
disp('Data Matrix:');
disp(data_matrix);
