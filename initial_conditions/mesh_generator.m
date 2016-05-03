function mesh_generator(f, xPoints, yPoints, fileName)

fid = fopen(fileName, 'w');

fprintf(fid, 'x, y, Temperature\n');

for x = xPoints
    for y = yPoints
        fprintf(fid, '%f, %f, %f\n', x, y, f(x, y));
    end
end

fclose(fid);

end