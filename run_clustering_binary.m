% Cluster the specimens into 2 groups based on our inconsistent surface mapping result
% 
% Reference:
% G. P. T. Choi, D. Qiu and L. M. Lui, 
% "Shape analysis via inconsistent surface registration."
% Proceedings of the Royal Society A, 476(2242), 20200147, 2020.
% 
% Copyright (c) 2020, Gary Pui-Tung Choi, Di Qiu, Lok Ming Lui

addpath('code')

%% phylogenetic groupings
[num, txt] = xlsread('HDM2016_genera.xlsx');
label_gt = zeros(50,1);
for i = 1:50
    if strcmp(txt{i,2},'Alouatta') 
        label_gt(i) = 1; %o
    elseif strcmp(txt{i,2},'Ateles') 
        label_gt(i) = 2; %s
    elseif strcmp(txt{i,2},'Brachyteles') 
        label_gt(i) = 3; %^
    elseif strcmp(txt{i,2},'Callicebus') 
        label_gt(i) = 4; %<
    elseif strcmp(txt{i,2},'Saimiri') 
        label_gt(i) = 5; %d
    end
end

%% binary clustering
load('inconsistent_registration_result.mat') % load the precomputed result

% insectivore vs non-insectivore
alpha = 0;
beta = 0;
gamma = 1;

% Folivore vs non-folivore
% alpha = 0.25;
% beta = 0.375;
% gamma = 0.375;

delta = alpha*mu_all + beta/2*Hndiff_all + gamma/2*Kndiff_all;
D = min(delta, delta');

%  hierarchical clustering
Z = linkage(D);
IDX = cluster(Z,'maxclust',2);

% visualization on the MDS plane 
C = mdscale(D,2);
colors = linspecer(max(IDX));
figure;
hold on;
for t = 1:length(D)
    if label_gt(t) == 1
        s = 'o';
    elseif label_gt(t) == 2
        s = 's';
    elseif label_gt(t) == 3
        s = '^';
    elseif label_gt(t) == 4
        s = '<';
    else
        s = 'd';
    end
    plot(C(t,1),C(t,2),['k',s],'MarkerFaceColor',colors(IDX(t),:),'LineWidth',2,'MarkerSize',13);
    text(C(t,1)+0.0015,C(t,2)+0.0015, num2str(t),'FontSize',20); 
end
set(gca,'linewidth',4,'FontSize',22);
xlabel('MDS Coordinate 1')
ylabel('MDS Coordinate 2')
