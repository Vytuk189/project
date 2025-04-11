% Open the file
filename = 'flowdata_cylinder2_nondimensional_CPP_N100_Re80.000000_iter100000_CFL05.txt';
fileID = fopen(filename, 'r');

% Read the first line (simulation time)
time = str2double(fgetl(fileID));

% Initialize an empty matrix to hold the data
data = [];

% Read the file line by line
line = fgetl(fileID);
while ischar(line)
    % Remove any leading or trailing white spaces
    line = strtrim(line);
    
    % Find all occurrences of [ ] and extract them
    data_points = regexp(line, '\[([^\]]+)\]', 'match');  % Find all substrings within brackets
    
    % Loop through all the found data points (which are inside [])
    for i = 1:length(data_points)
        % Remove the brackets and split the string by spaces
        point_str = data_points{i};
        point_values = str2num(point_str(2:end-1)); % Convert string to numbers, excluding the brackets
        
        % Append this data point to the data matrix
        data = [data; point_values];
    end
    
    % Read the next line
    line = fgetl(fileID);
end

% Close the file
fclose(fileID);


% Reshape the data to extract x, y coordinates and pressure values
x = data(:, 4);  % X coordinates
y = data(:, 5);  % Y coordinates
pressure = data(:, 1);  % Pressure values
horizontal_velocity = data(:, 2);  % Horizontal velocity (second column)
vertical_velocity = data(:, 3);    % Vertical velocity (third column)
total_velocity = sqrt(horizontal_velocity.^2 + vertical_velocity.^2);

% Define grid boundaries (adjust the range of x and y to fit the simulation area)
x_min = min(x); x_max = max(x);
y_min = min(y); y_max = max(y);

% Create a grid of x and y values for the contour plot
[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));

% Interpolate the pressure data onto the grid
pressure_grid = griddata(x, y, pressure, x_grid, y_grid, 'cubic');

% % Create the contour plot
% figure;
% contourf(x_grid, y_grid, pressure_grid, 20, 'LineColor', 'none');  % 20 levels for contour
% colorbar;  % Display the color scale
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Pressure Contour Map');

% Interpolate the horizontal velocity onto the grid
horizontal_velocity_grid = griddata(x, y, horizontal_velocity, x_grid, y_grid, 'cubic');

% Interpolate the vertical velocity onto the grid
vertical_velocity_grid = griddata(x, y, vertical_velocity, x_grid, y_grid, 'cubic');

% Interpolate the total velocity onto the grid
total_velocity_grid = griddata(x, y, total_velocity, x_grid, y_grid, 'cubic');

% % Plot Horizontal Velocity Contour Map
% figure;
% contourf(x_grid, y_grid, horizontal_velocity_grid, 20, 'LineColor', 'none');
% colorbar;
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Horizontal Velocity Contour Map');
% 
% % Plot Vertical Velocity Contour Map
% figure;
% contourf(x_grid, y_grid, vertical_velocity_grid, 20, 'LineColor', 'none');
% colorbar;
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Vertical Velocity Contour Map');
% 
% % Plot Total Velocity (Speed) Contour Map
% figure;
% contourf(x_grid, y_grid, total_velocity_grid, 20, 'LineColor', 'none');
% colorbar;
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Total Velocity Contour Map');


% Create a grid of x and y values for the streamline plot
[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));

% Interpolate the horizontal and vertical velocities onto the grid
horizontal_velocity_grid = griddata(x, y, horizontal_velocity, x_grid, y_grid, 'cubic');
vertical_velocity_grid = griddata(x, y, vertical_velocity, x_grid, y_grid, 'cubic');

% Create a figure for streamline plot
figure;
hold on;

% Plot streamlines using the interpolate velocity fields
n = 500;
startx = zeros(n,1);  % Seed points for streamlines in x-direction
starty = linspace(y_min, y_max, n);  % Seed points for streamlines in y-direction

% Create streamlines based on seed points
streamline(x_grid, y_grid, horizontal_velocity_grid, vertical_velocity_grid, startx, starty);

% Customize the plot
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Streamline Plot of Fluid Flow');
axis tight;
hold off;
