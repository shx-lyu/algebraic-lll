% main function: comparing the performance of algebraic LLL (bases in the complex field) and
% non-algebraic (bases in the real field) LLL

% Setting: bases are random integers; algebraic LLL algorithms are defined over rings: Z[\omega],Z[i],Z[\sqrt{-2}]
% author: Shanxiang Lyu, shanxianglyu@gmail.com, Jinan University

clc;tic
clear all;close all; warning off;
linestyles = cellstr(char('-','--','-','--','-','--'));
SetColors=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0    0.4470    0.7410
    0.8500    0.3250    0.0980];
Markers=['o','o','x','x','+','+'];
FIT=[];TIME=[];

Dimensions=[2:2:12];
Algorithms=[1:6];
 
legendbox={'$Z[\omega]$-ALLL','$Z[\omega]$-RLLL',...
   '$Z[i]$-ALLL','$Z[i]$-RLLL',...
   '$Z[\sqrt{-2}]$-ALLL','$Z[\sqrt{-2}]$-RLLL'};
 
 for n=Dimensions 
        for monte=1:1e2
  
            B=randi([0,4*n],n,n)+sqrt(-1)*randi([0,4*n],n,n);
            for j=Algorithms
                switch j
                case 1 
                expression = '[BR]=ALLL(B,3);';
                case 2 
                expression = '[BR]=RLLL(B,3);';
                case 3 
                expression = '[BR]=ALLL(B,1);';
                case 4 
                expression = '[BR]=RLLL(B,1);';
                case 5 
                expression = '[BR]=ALLL(B,2);';
                case 6 
                expression = '[BR]=RLLL(B,2);';
                end

            tic
            eval(expression);
            time_j(monte,j)=toc;
            fit_j(monte,j)= norm(BR(:,1));

            end
        end
        
        FIT=[FIT,mean(fit_j,1)'];
        TIME=[TIME,mean(time_j,1)'];

end

figure(1)
    for j=Algorithms
semilogy(Dimensions,FIT(j,:),[linestyles{j} Markers(j)],'Color',SetColors(j,:),'Linewidth',2);
        hold on;
        grid on;
    end
hold off;
legend(legendbox(Algorithms),'Interpreter','latex');
xlabel('Dimension n');ylabel('Length of the Shortest Vector');
 
figure(2)
    for j=Algorithms
semilogy(Dimensions,TIME(j,:),[linestyles{j} Markers(j)],'Color',SetColors(j,:),'Linewidth',2);
        hold on;
        grid on;
    end
hold off;
legend(legendbox(Algorithms),'Interpreter','latex');
xlabel('Dimension n');ylabel('Running Time/s');
