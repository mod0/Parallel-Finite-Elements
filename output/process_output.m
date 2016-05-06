N=100;
M = 100;


data = csvread('data.csv.10');
data = sortrows(data, [1 2]);
surf(reshape(data(:,1),N,M), reshape(data(:,2),N,M), reshape(data(:,3),N,M))
xlabel('x');
ylabel('y');
zlabel('Temperature');