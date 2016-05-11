times2 = [36.905 18.405];
times4 = [62.973 19.962];
times8 = [72.457 10.927];

y = [times2; times4; times8];

threads = [2 4 6];

bar(threads, y);
set(gca,'XTickLabel',{'2', '4', '8'});
xlabel('Number of Subdomains/Threads');
ylabel('Times (s)');

legend('Serial', 'Parallel');
title('Running Times on 128x128 Domain');