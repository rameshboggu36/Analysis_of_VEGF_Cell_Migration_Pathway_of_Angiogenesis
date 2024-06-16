close all;
clear all;
clc;
y0 = [100 80 70 60 50 0 0 40 0 30 0 20 0 10 0 50 0 0 20 0 40 0 10 0];
tspan = [0 40];
[t,y] = ode45(@vegf_up,tspan,y0);

D = {'VEGF', 'VEGF-R', 'VR', '(VEGF-VEGFR)2', 'NCK', 'NCK_A', 'NCKA_VR2', 'PAK2', 'PAK2_A', 'P38',...
    'P38_A', 'MAPKAPK2', 'MAPKAPK2_A', 'HSP27', 'HSP27_A', 'SHB', 'SHB_A', 'SHB_A_VR2',...
    'PRAK', 'PRAK_A', 'FAK', 'FAK_A', 'PAXILLIN', 'PAXILLIN_A'};

disp('Variation of Concentrations of all proteins Over the time')
% [1,2,3,4,5,6,16,17]
for i = 1:24
figure(i);
plot(t,y(:,i),'linewidth', 2)
ylabel(D(i))
xlabel('Time(sec)')
str = "Dynamics of "+D(i);
title(str)
file = D(i)+".jpg";
saveas(figure(i),file ,'jpg');
end 

disp('PCA using the inbuilt method using protein dynamics')
featureNames = D(1:24);

[coeff, ~, ~, ~, ~] = pca(y);
for i = 1:size(coeff, 2)
    [~, idx] = sort(coeff(:, i), 'descend');
end

% Display the sum of loadings for each feature along with feature names, sorted in ascending order
totalLoadings = sum(coeff, 2);
[~, sortedIndices] = sort(abs(totalLoadings));
fprintf('\nSum of Loadings for Each Feature:\n');
for i = 1:length(sortedIndices)
    fprintf('%s %f\n', featureNames{sortedIndices(i)}, totalLoadings(sortedIndices(i)));
end

