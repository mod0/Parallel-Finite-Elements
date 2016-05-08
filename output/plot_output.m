function  plot_output(data)
sortrows(data, [1 2]);

x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

m = length(unique(y));
n = length(unique(x));

surf(reshape(x, m, n), reshape(y, m, n), reshape(z, m, n));
xlabel('x');
ylabel('y');
zlabel('Temperature');

end