function time_mesh_generator(f, tPoints, xPoints, yPoints, fileName)

for i = 1:length(tPoints)
    fxy = @(x, y) f(tPoints(i), x, y);
    mesh_generator(fxy, xPoints, yPoints, [fileName '.' num2str(i-1)]);
end

end