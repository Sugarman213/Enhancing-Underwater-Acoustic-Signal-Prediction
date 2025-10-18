%% 优化VMD-CNN-BiLSTM-Attention

disp(' ')
disp('')

%% 

clc;
clear 
close all
addpath(genpath(pwd))


load VMD_IMF_SHIP1.mat
X=Co_IMF1';


num_samples = length(X);                         
kim = 10;                      
zim =  1;                     
or_dim = size(X,2);

%  
for i = 1: num_samples - kim - zim + 1
    res(i, :) = [reshape(X(i: i + kim - 1,:), 1, kim*or_dim), X(i + kim + zim - 1,end)];
end


% 
outdim = 1;                                  
num_size = 0.8;                              
num_train_s = round(num_size * num_samples); 
f_ = size(res, 2) - outdim;                


P_train = res(1: num_train_s, 1: f_)';
T_train = res(1: num_train_s, f_ + 1: end)';
M = size(P_train, 2);

P_test = res(num_train_s + 1: end, 1: f_)';
T_test = res(num_train_s + 1: end, f_ + 1: end)';
N = size(P_test, 2);

%  
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

%%  数据平铺 %% 
for i = 1:size(P_train,2)
    trainD{i,:} = (reshape(p_train(:,i),size(p_train,1),1,1));
end

for i = 1:size(p_test,2)
    testD{i,:} = (reshape(p_test(:,i),size(p_test,1),1,1));
end


targetD =  t_train;
targetD_test  =  t_test;

numFeatures = size(p_train,1);


%% 
popsize=10;  
maxgen=20;   
fobj = @(x)objective(x,f_,X,ps_output);
% 优化参数设置
lb = [1 200];
ub = [20 5000];    
dim = length(lb);

% 'DBO','GWO','OOA','PSO','SABO','SCSO','SSA','BWO','RIME','WOA','HHO','NGO';

[Best_score,Best_pos,curve]=BWO(popsize,maxgen,lb,ub,dim,fobj); 
setdemorandstream(pi);

%% 
figure
plot(curve,'r-','linewidth',2)
xlabel('进化代数')
ylabel('均方误差')
legend('最佳适应度')
title('进化曲线')

%% 
[~,optimize_T_sim] = objective(Best_pos,f_,X,ps_output);
setdemorandstream(pi);
%% 
str={'真实值','优化VMD-CNN-BiLSTM-Attention'};
figure('Units', 'pixels', ...
    'Position', [300 300 860 370]);
plot(T_test,'-','Color',[0.8500 0.3250 0.0980]) 
hold on
plot(optimize_T_sim,'-.','Color',[0.4940 0.1840 0.5560]) 
legend(str)
set (gca,"FontSize",12,'LineWidth',1.2)
box off
legend Box off


%% 
test_y = T_test;
Test_all = [];

y_test_predict = optimize_T_sim;
[test_MAE,test_MAPE,test_MSE,test_RMSE,test_R2]=calc_error(y_test_predict,test_y);
Test_all=[Test_all;test_MAE test_MAPE test_MSE test_RMSE test_R2];
save SHIPCO1_OP_VMD_CNN_BiLSTM_Attention optimize_T_sim



% str1=str(2:end);
% str2={'MAE','MAPE','MSE','RMSE','R2'};
% data_out=array2table(Test_all);
% data_out.Properties.VariableNames=str2;
% data_out.Properties.RowNames=str1;
% disp(data_out)
