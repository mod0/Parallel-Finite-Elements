f = @(t, x, y) sin(t + x + y^2);

t = 1:5;
x = 1:10;
y = 1:10;

time_mesh_generator(f, t, x, y, 'test.csv');