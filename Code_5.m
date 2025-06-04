
clear

rng('default')      % We save the seed
state = rng;


Code_1_German
%Code_1_Student



Num = 5;                               % Number of Folds.
M=floor(N/Num);                          % Approximate number of data in test set.
VecL = logspace(-3,8,60);              % We define the grid of lambda.
F = 1; % If F=1 then the fairness constraint is enforced in the train set and otherwise in the test set.

Name = sprintf('Results_May24_GERMAN_Test_aux');

indexChange = randperm(N);

y0=y0(indexChange);
X0=X0(indexChange,:);
isnotSens=isnotSens(indexChange);

isSens=isSens(indexChange);

ySens=y0(isSens);
yNoSens=y0(~isSens);

t0_Tot=tic;

a_ini=1e-2;
b_ini=1e-2;


save('Results_May24_GERMAN_Test_Global_aux.mat')

for n=1:Num

    [X,y,z,XT,yT,zT,Ind_Train,Ind_Test,medSensF,medNotSensF,medSensFT,medNotSensFT]=Code_2_fun_CV(X0,y0,isnotSens,M,n);
    
    iter = sprintf('_CV_%d.mat', n);

    if det(X0'*X0)==0 || det(X'*X)==0 
        disp('Error: Null determinant')
        break
    end
    
    t0_iter=tic;
    Code_4_fun(X,y,z,XT,yT,zT,F,VecL,strcat(Name,iter))
    t1_iter(n)=toc(t0_iter)

end

t1_Tot=toc(t0_Tot)

save('Results_May24_GERMAN_Test_CV_Full_aux.mat')











