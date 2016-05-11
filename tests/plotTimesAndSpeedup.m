function plotTimesAndSpeedup(threads, serialTime, parallelTimes, domainSize)
    figure; hold on;
    plot(threads, parallelTimes, 'O-');
    plot(threads, serialTime * ones(size(threads)));
    xlabel('Number of Threads');
    ylabel('Time (s)');
    title(sprintf('Running Time vs Number of Threads on %dx%d Domain', domainSize, domainSize));
    legend('Parallel Times', 'Serial Time');
    
    figure; hold on;
    plot(threads, serialTime ./ parallelTimes, 'O-');
    xlabel('Number of Threads');
    ylabel('Speedup');
    title(sprintf('Speedup on %dx%d Domain', domainSize, domainSize));
end