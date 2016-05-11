numThreads = [1 2 4 8 16 32];

serial64 = 0.284;
parallel64 = [0.238, 0.663, 0.736, 0.448, 0.235, 0.218];
plotTimesAndSpeedup(numThreads, serial64, parallel64, 64);

seria1128 = 2.544;
parallel128 = [2.529 15.122 10.415 6.055 3.005 2.121];
plotTimesAndSpeedup(numThreads, seria1128, parallel128, 128);

seria1256 = 39.321;
parallel256 = [38.851 150.077 162.433 90.035 41.760 27.062];
plotTimesAndSpeedup(numThreads, seria1256, parallel256, 256);

seria1512 = 598.103;
parallel512 = [601.261 2335.712 2402.341 1340.412 620.412 332.220];
plotTimesAndSpeedup(numThreads, seria1512, parallel512, 512);