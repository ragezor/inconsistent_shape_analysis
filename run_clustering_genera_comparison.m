% Comparison: classification into five groups using traditional methods
% 
% Reference:
% G. P. T. Choi, D. Qiu and L. M. Lui, 
% "Shape analysis via inconsistent surface registration."
% Proceedings of the Royal Society A, 476(2242), 20200147, 2020.
% 
% Copyright (c) 2020, Gary Pui-Tung Choi, Di Qiu, Lok Ming Lui

addpath('code')
addpath('code/natsortfiles')
addpath('code/Comparison')

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

%% Procrustes and ICP
prefix = 'HDM 2016_preprocessed/';

data = dir([prefix,'*.mat']);
data_list = natsortfiles({data.name})';

distance_procrustes = zeros(length(data_list));
distance_icp = zeros(length(data_list));
for i = 1:length(data_list)
    %% load source data and preprocess
    disp(['..................i = ',num2str(i)]);
    
    source_mat =  load([prefix,data_list{i}]);
    source_face = source_mat.f;
    source_vertex = source_mat.v;
    source_gauss_curvature = source_mat.K;
    source_mean_curvature = source_mat.H;
    landmark_source_index = source_mat.lm;

    %% Target
    for j = 1:length(data_list)
        %% load target data and preprocess
        target_mat = load([prefix,data_list{j}]);

        target_face = target_mat.f;
        target_vertex = target_mat.v;
        target_gauss_curvature = target_mat.K;
        target_mean_curvature = target_mat.H;
        landmark_target_index = target_mat.lm;
        
        [distance_procrustes(i,j)] = procrustes(target_vertex(landmark_target_index,:),source_vertex(landmark_source_index,:));
        [~, ~, ER] = icp_matlabcentral(target_vertex',source_vertex',10);
        distance_icp(i,j) = ER(end);
    end
end

%% clustering into 5 groups
% load the precomputed result
load('dissimilarity_comparison.mat');

D = (distance_procrustes+distance_procrustes')/2;
% D = (distance_icp+distance_icp')/2;
D = D-diag(diag(D));

% k-means
C = mdscale(D,2); 
IDX = kmeans(C,5,'Replicates',100);

count = 0;
for i = 1:length(D)
    for j = i+1:length(D)
        if (label_gt(i) == label_gt(j) && IDX(i) == IDX(j)) || (label_gt(i) ~= label_gt(j) && IDX(i) ~= IDX(j))
            % correct
            count = count + 1;
        end
    end
end
accuracy = count/((length(D)-1)*(length(D))/2);

% visualization on the MDS plane 
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
    text(C(t,1)+0.001,C(t,2)+0.001, num2str(t),'FontSize',20); % for Procrustes
%     text(C(t,1)+0.000025,C(t,2)+0.000025, num2str(t),'FontSize',20); % for ICP
end
set(gca,'linewidth',4,'FontSize',22);
xlabel('MDS Coordinate 1')
ylabel('MDS Coordinate 2')
