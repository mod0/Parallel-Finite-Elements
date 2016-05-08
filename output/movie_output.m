function movieFrames = movie_output(data)

catData = vertcat(data{:});
z = catData(:, 3);
zMin = min(z);
zMax = max(z);

numSteps = length(data);
movieFrames(numSteps) = struct('cdata', [], 'colormap', []);

figure;

for i = 1:numSteps
    plot_output(data{i});
    zlim([zMin, zMax]);
    drawnow;
    movieFrames(i) = getframe;
end

end