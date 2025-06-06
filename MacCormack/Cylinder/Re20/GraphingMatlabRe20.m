Re = 10;
N = 100;

% Open the file
filename = 'flowdata_cylinder_NR25.000000_Re20.000000_Iter799000.txt';
fileID = fopen(filename, 'r');

% % Read the first line (simulation time)
% time = str2double(fgetl(fileID));

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
structure = data(:, 6);
if structure == 0
    horizontal_velocity = NaN;  % Horizontal velocity (second column)
    vertical_velocity = NaN;
end


% Define grid boundaries (adjust the range of x and y to fit the simulation area)
x_min = min(x); x_max = max(x);
y_min = min(y); y_max = max(y);

% Create a grid of x and y values for the contour plot
[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 200), linspace(y_min, y_max, 200));


% You may need to adjust the dimensions of your grid (100x100) depending on the size of your data
structure_grid = griddata(x, y, structure, x_grid, y_grid);









% Interpolate the pressure data onto the grid
pressure_grid = griddata(x, y, pressure, x_grid, y_grid, 'cubic');

% Create the contour plot
figure;
contourf(x_grid, y_grid, pressure_grid, 20, 'LineColor', 'none');  % 20 levels for contour
% Create a custom colormap from dark blue to dark red
colorbar;  % Display the color scale
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Pressure Contour Map');

% Interpolate the horizontal velocity onto the grid
horizontal_velocity_grid = griddata(x, y, horizontal_velocity, x_grid, y_grid, 'cubic');

% Interpolate the vertical velocity onto the grid
vertical_velocity_grid = griddata(x, y, vertical_velocity, x_grid, y_grid, 'cubic');

% Interpolate the total velocity onto the grid
total_velocity_grid = griddata(x, y, total_velocity, x_grid, y_grid, 'cubic');

% Plot Horizontal Velocity Contour Map
figure;
contourf(x_grid, y_grid, horizontal_velocity_grid, 20, 'LineColor', 'none');
colorbar;
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Horizontal Velocity Contour Map');

% Plot Vertical Velocity Contour Map
figure;
contourf(x_grid, y_grid, vertical_velocity_grid, 20, 'LineColor', 'none');
colorbar;
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Vertical Velocity Contour Map');

% Plot Total Velocity (Speed) Contour Map
figure;
contourf(x_grid, y_grid, total_velocity_grid, 20, 'LineColor', 'none');
colorbar;
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Total Velocity Contour Map');


% Create a grid of x and y values for the streamline plot
[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100));

% Interpolate the horizontal and vertical velocities onto the grid
horizontal_velocity_grid = griddata(x, y, horizontal_velocity, x_grid, y_grid, 'cubic');
vertical_velocity_grid = griddata(x, y, vertical_velocity, x_grid, y_grid, 'cubic');


% Plot streamlines using the interpolate velocity fields
n = 40;
startx = zeros(1,n);  % Seed points for streamlines in x-direction
starty = linspace(y_min, y_max, n);  % Seed points for streamlines in y-direction


% Cylinder parameters
cx = 0.8;   % X-coordinate of the cylinder center
cy = 0.5;   % Y-coordinate of the cylinder center
r = 0.2;   % Radius of the cylinder

% Define the angle range for the right half of the cylinder (edge)
theta = linspace(-4*pi/10, 4*pi/10, 10);  % Generate 20 points along the right edge
% Calculate the coordinates of the seed points on the edge of the cylinder
seedX = cx + 1.02*r * cos(theta);  % x-coordinates of the seed points
seedY = cy + 1.02*r * sin(theta);  % y-coordinates of the seed points
% Add seed points from the cylinder's right edge to the startx and starty arrays
startx = [startx, seedX];  % Concatenate seedX to startx (make sure to transpose seedX if needed)
starty = [starty, seedY];  % Concatenate seedY to starty (make sure to transpose seedY if needed)

m = 10;
vortexx = linspace(cx, cx+1*r, m);
vortexy = ones(1,m)*(cy+0.1);

startx = [startx, vortexx];  % Concatenate seedX to startx (make sure to transpose seedX if needed)
starty = [starty, vortexy];  % Concatenate seedY to starty (make sure to transpose seedY if needed)

vortexx = linspace(cx, cx+1*r, m);
vortexy = ones(1,m)*(cy-0.1);

startx = [startx, vortexx];  % Concatenate seedX to startx (make sure to transpose seedX if needed)
starty = [starty, vortexy];  % Concatenate seedY to starty (make sure to transpose seedY if needed)

k = 5;
edgex = linspace(1.05, 2, k);
edgey = zeros(1,k)+0.001;
startx = [startx, edgex];  % Concatenate seedX to startx (make sure to transpose seedX if needed)
starty = [starty, edgey];  % Concatenate seedY to starty (make sure to transpose seedY if needed)
edgex = linspace(1.05, 2, k);
edgey = 1*ones(1,k)-0.001;
startx = [startx, edgex];  % Concatenate seedX to startx (make sure to transpose seedX if needed)
starty = [starty, edgey];  % Concatenate seedY to starty (make sure to transpose seedY if needed)


%WITH SQUARES

%Create streamlines based on seed points

% Create a new figure
figure;

% Set the figure background color to white
set(gcf, 'Color', 'w');  % Set the figure background color to white

% Plot the structure data
imagesc(x_grid(1,:), y_grid(:,1), structure_grid);
colormap([0 0 0; 1 1 1]); % Black for 0, White for 1
axis xy;  % Ensure the y-axis is in the correct orientation
hold on;  % Keep the current plot and allow subsequent plots

h = streamline(x_grid, y_grid, horizontal_velocity_grid, vertical_velocity_grid, startx, starty);

% Set streamline properties (make them black)
set(h, 'Color', 'k');  % 'k' stands for black



% Set axis to 1:1 ratio (no distortion of shapes)
axis equal;  % Ensure scaling is 1:1



% Tighten the axes
axis tight;

xlim([0 2]);
ylim([0 1]);

% Customize the plot
xlabel('X', 'Color', 'k', 'FontSize', 12,'FontWeight', 'bold');  % Set the x-axis label color to black
ylabel('Y', 'Color', 'k', 'FontSize', 12,'FontWeight', 'bold');  % Set the y-axis label color to black

% Customize the title
title(['Proudnice pro N = ' num2str(N) ' a Re = ' num2str(Re)], ...
    'Color', 'k', ...             % Set title color to black
    'FontSize', 16, ...           % Make the title a bit bigger
    'VerticalAlignment', 'bottom', ... % Move the title a bit further away
    'FontWeight', 'bold');        % Make the title bold

% Set the axis ticks to be black (this is usually default, but just to ensure)
set(gca, 'XColor', 'k', 'YColor', 'k');  % Set both X and Y axis ticks to black


hold off;

% Define the file name
filename = ['StreamlineN' num2str(N) 'Re' num2str(Re) '.png'];

% Save the figure to a file automatically as PNG
saveas(gcf, filename, 'png');


% % WITH DOTS
% % Create a new figure
% figure;
% 
% % Plot the structure data as black and white dots
% % Find indices of points where structure is 0 and 1
% indices_zeros = structure_grid == 0;
% indices_ones = structure_grid == 1;
% 
% % Scatter plot for structure 0 (black dots)
% scatter(x_grid(indices_zeros), y_grid(indices_zeros), 10, 'k', 'filled'); % Black dots for 0
% hold on;
% 
% h = streamline(x_grid, y_grid, horizontal_velocity_grid, vertical_velocity_grid, startx, starty);
% 
% % Set streamline properties (make them black)
% set(h, 'Color', 'k');  % 'k' stands for black
% 
% % Customize the plot
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Streamline Plot of Fluid Flow');
% 
% % Set axis to 1:1 ratio (no distortion of shapes)
% axis equal;  % Ensure scaling is 1:1
% 
% % Set background to white
% set(gca, 'Color', 'w');  % Set the background color to white
% 
% % Tighten the axes
% axis tight;
% 
% % Hold off to finish the plotting
% hold off;
