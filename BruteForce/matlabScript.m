clc;
close all;
clearvars;
format long;

fstring = './BruteForce';
probName = 'BruteForce';

runs = 3;
procs = 5;
prob = 13;

problem_size = [1 2 4 8 16 32 64 128 256 512 1024 2048 4096];
FLOPs = zeros(size(problem_size));
for i = 1:length(problem_size)
    FLOPs(i) = flops(problem_size(i));
end
processor_size = [0 2 4 8 16];

alg_time = zeros(length(processor_size),length(problem_size));
e2e_time = zeros(length(processor_size),length(problem_size));

W = readmatrix('data.txt');
display(W);

for i = 1:length(W)
    prob_ind = find(problem_size==W(i,1));
    proc_ind = find(processor_size==W(i,2));
    display(['N = ',num2str(prob_ind),' L = ',num2str(proc_ind)]);
    alg_time(proc_ind,prob_ind) = alg_time(proc_ind,prob_ind) + W(i,4);
    e2e_time(proc_ind,prob_ind) = e2e_time(proc_ind,prob_ind) + W(i,5);
end

alg_time = alg_time./runs;

e2e_time = e2e_time./runs;

serial_time = alg_time(1,:);

speedup = alg_time./serial_time;
speedup = 1./speedup;

display(speedup);
processor_size(1) = 1;
efficiency = speedup'./processor_size;
efficiency = efficiency';

legendArray = "Threads = "+string(processor_size(1,2:end));
legendArray = ["Serial",legendArray];

figure(1);
for i = 1:length(processor_size)
    semilogy(log2(problem_size),alg_time(i,:),'o-','LineWidth',2,'MarkerSize',2.5,'MarkerFaceColor','auto');
    hold on;
end
grid on;
title(['Algorithm Run Time - ',probName]);
legend(legendArray,'Location','best');
xlabel('logN');
ylabel('Run Time (in seconds)');
saveas(gca,[fstring,'_AlgTime'],'png');

figure(2);
for i = 1:length(processor_size)
    plot(log2(problem_size),speedup(i,:),'o-','LineWidth',2,'MarkerSize',2.5,'MarkerFaceColor','auto');
    hold on;
end
grid on;
legend(legendArray,'Location','best');
title(['Speed Up - ',probName]);
xlabel('logN');
ylabel('Speed Up');
saveas(gca,[fstring,'_SpeedUp_Prob'],'png');

figure(3);
for i = 8:length(problem_size)
    plot(processor_size,speedup(:,i),'o-','LineWidth',2,'MarkerSize',2.5,'MarkerFaceColor','auto');
    hold on;
end
grid on;
legend("N = "+string(problem_size(1,8:end)),'Location','northwest');
title(['Speed Up - ',probName]);
ylabel('Speed Up');
xlabel('Threads');
saveas(gca,[fstring,'_SpeedUp_Proc'],'png');

figure(4);
for i = 1:length(processor_size)
    plot(log2(problem_size),efficiency(i,:),'o-','LineWidth',2,'MarkerSize',2.5,'MarkerFaceColor','auto');
    hold on;
end
grid on;
legend(legendArray,'Location','northwest');
title(['Efficiency - ',probName]);
ylabel('Effiecncy (Speed Up / Threads)');
xlabel('logN');
saveas(gca,[fstring,'_Effi'],'png');

figure(5);
for i = 1:length(processor_size)
    plot(log2(problem_size),FLOPs./(alg_time(i,:)*10^9),'o-','LineWidth',2,'MarkerSize',2.5,'MarkerFaceColor','auto');
    hold on;
end
grid on;
title(['Performance - ',probName]);
legend(legendArray,'Location','best');
xlabel('logN');
ylabel('Performance GFLOPS');
saveas(gca,[fstring,'_Perf'],'png');

function f = flops(n)
    if n <= 32
        f = 2*n*n*n;
    else
        f = 7*flops(n/2) + 15*n*n/4;
    end
end
