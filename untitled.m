

% --- setup grid ---
nx = 200;  % number of points in x-direction
ny = 150;  % number of points in y-direction
x = linspace(-5, 5, nx);  % nx points equally spaced between -5...5
y = linspace(-6, 6, ny);  % ny points equally spaced between -6...6
[X, Y] = ndgrid(x, y);  % 2D array (matrix) of points across x and y
Z = zeros(nx, ny);  % initialize output of size (nx, ny)

% --- evaluate across grid ---
for i = 1:nx
    for j = 1:ny
        Z(i, j) = func(X(i, j), Y(i, j));
    end
end

% --- contour plot ---
figure();  % start a new figure
contour(X, Y, Z, 50);  % using 50 contour lines.  
colorbar();  % add a colorbar
xlabel('x');  % labels for axes
ylabel('y');

function [z] = func(x, y)
z = 3*(x-2)^2 + (y+1)^2;
end