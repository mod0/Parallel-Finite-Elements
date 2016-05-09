data = process_output();
movieFrames = movie_output(data);

movie(movieFrames, 100, 60)

writerObj = VideoWriter('movie.avi');
open(writerObj);
writeVideo(writerObj, movieFrames);
close(writerObj);