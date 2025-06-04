
clear
Num=5;

for n=1:Num
    clearvars -except Num n Tasa_Exito_CV_Train_aux Tasa_Exito_CV_Test_aux Dif_CV_Train_aux Dif_CV_Test_aux M_p97_5_TE_Train M_p97_5_TE_Test M_p2_5_TE_Train M_p2_5_TE_Test M_H_beta_Train M_H_beta_Test
    filename = sprintf('Results_May24_GERMAN_Train_1_CV_%d.mat', n);
    load(filename);
    

    for i = 1:length(VecL)
        beta = M_beta{i};
        Sigma = M_Sigma{i};

        Nsim=10^4;
        N = length(y);
        NT = length(yT);

        beta_extr = mvnrnd(beta,Sigma,Nsim)';
        py=1./(1+exp(-beta_extr'*X'));

        H_beta_Train_aux = beta_extr'*z*z'*beta_extr;
        H_beta_Train(:,i) = diag(H_beta_Train_aux);

        H_beta_Test_aux = beta_extr'*zT*zT'*beta_extr;
        H_beta_Test(:,i) = diag(H_beta_Test_aux);


        % % M_Tasa_Exito_Train_aux=zeros(N,Nsim);
        % % 
        % % for j=1:Nsim
        % %     for k=1:N
        % %         if py(j,k)>=0.5
        % %             if y(k)==1
        % %                 M_Tasa_Exito_Train_aux(k,j)=1;
        % %             end
        % %         else
        % %             if y(k)==-1
        % %                 M_Tasa_Exito_Train_aux(k,j)=1;
        % %             end
        % %         end
        % %     end
        % % 
        % %     V_Tasa_Exito_Train(j)=mean(M_Tasa_Exito_Train_aux(:,j));
        % % end
        % % 
        % % 
        % % M_Tasa_Exito_Train{i}=V_Tasa_Exito_Train;
        % % 
        % % Tasa_Exito_Train(i) = mean(V_Tasa_Exito_Train);
        % % p97_5_Train(i)=prctile(V_Tasa_Exito_Train,97.5);
        % % p2_5_Train(i)=prctile(V_Tasa_Exito_Train,2.5);
        % % 
        % % Dif_Train(i) = 2*trace(z'*Sigma*z)+beta'*z*z'*beta;
        % % 
        % % py=1./(1+exp(-beta_extr'*XT'));
        % % M_Tasa_Exito_Test_aux=zeros(NT,Nsim);
        % % for j=1:Nsim
        % %     for k=1:NT
        % %         if py(j,k)>=0.5
        % %             if yT(k)==1
        % %                 M_Tasa_Exito_Test_aux(k,j)=1;
        % %             end
        % %         else
        % %             if yT(k)==-1
        % %                 M_Tasa_Exito_Test_aux(k,j)=1;
        % %             end
        % %         end
        % %     end
        % % 
        % %     V_Tasa_Exito_Test(j)=mean(M_Tasa_Exito_Test_aux(:,j));
        % % end
        % % M_Tasa_Exito_Test{i}=V_Tasa_Exito_Test;
        % % 
        % % Tasa_Exito_Test(i) = mean(V_Tasa_Exito_Test);
        % % p97_5_Test(i)=prctile(V_Tasa_Exito_Test,97.5);
        % % p2_5_Test(i)=prctile(V_Tasa_Exito_Test,2.5);
        % % 
        % % Dif_Test(i) = 2*trace(zT'*Sigma*zT)+beta'*zT*zT'*beta;

    end


    M_H_beta_Train{n}=H_beta_Train;
    M_H_beta_Test{n}=H_beta_Test;


    M_p97_5_TE_Train(:,n)=p97_5_Train';
    M_p2_5_TE_Train(:,n)=p2_5_Train';
    M_p97_5_TE_Test(:,n)=p97_5_Test';
    M_p2_5_TE_Test(:,n)=p2_5_Test';



    Tasa_Exito_CV_Train_aux(:,n)=Tasa_Exito_Train';
    Tasa_Exito_CV_Test_aux(:,n)=Tasa_Exito_Test';
    Dif_CV_Train_aux(:,n)=Dif_Train';
    Dif_CV_Test_aux(:,n)=Dif_Test';

end

for i=1:length(VecL)
    Tasa_Exito_CV_Train(i) = mean(Tasa_Exito_CV_Train_aux(i,:));
    Tasa_Exito_CV_Test(i) = mean(Tasa_Exito_CV_Test_aux(i,:));
    Dif_CV_Train(i) = mean(Dif_CV_Train_aux(i,:));
    Dif_CV_Test(i) = mean(Dif_CV_Test_aux(i,:));

    p97_5_Train(i) = mean(M_p97_5_TE_Train(i,:));
    p2_5_Train(i) = mean(M_p2_5_TE_Train(i,:));
    p97_5_Test(i) = mean(M_p97_5_TE_Test(i,:));
    p2_5_Test(i) = mean(M_p2_5_TE_Test(i,:));

end

%save('Results_Abr24_GERMAN_Test_MA_NoRidge_Complete_3.mat')

%% Gráficas

%load('Results_May24_STUDENT_Test_Complete_Final_1.mat')
load('Results_May24_GERMAN_Test_Complete_Final_1.mat')

Lim = length(VecL);

figure (1)
semilogx(VecL(1:Lim),Tasa_Exito_CV_Train(1:Lim),'r','linewidth',2)
hold on
semilogx(VecL(1:Lim),Tasa_Exito_CV_Test(1:Lim),'b','linewidth',2)
semilogx(VecL(1:Lim),p97_5_Train(1:Lim),'--r','linewidth',2)
semilogx(VecL(1:Lim),p2_5_Train(1:Lim),'--r','linewidth',2)
semilogx(VecL(1:Lim),p97_5_Test(1:Lim),'--b','linewidth',2)
semilogx(VecL(1:Lim),p2_5_Test(1:Lim),'--b','linewidth',2)
set(gca,'FontSize',14)
xlabel('$\lambda$','Interpreter','latex','FontSize',24)
ylabel('Success Rate','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','95\% PI (Train)','','95\% PI (Test)','Interpreter','latex','FontSize',20)
grid on
grid minor


figure (2)
semilogx(VecL(1:Lim),Tasa_Exito_CV_Train(1:Lim),'r','linewidth',2)
hold on
semilogx(VecL(1:Lim),Tasa_Exito_CV_Test(1:Lim),'b','linewidth',2)
set(gca,'FontSize',14)
xlabel('$\lambda$','Interpreter','latex','FontSize',24)
ylabel('Success Rate','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','Interpreter','latex','FontSize',20)
grid on
grid minor


figure (3)
loglog(VecL(1:Lim),Dif_CV_Train(1:Lim),'--r','linewidth',2)
hold on
loglog(VecL(1:Lim),Dif_CV_Test(1:Lim),'--b','linewidth',2)
grid on
grid minor
set(gca,'FontSize',14)
set(gca,'FontSize',14)
xlabel('$\lambda$','Interpreter','latex','FontSize',24)
ylabel('$\mathcal{F}_0$','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','Interpreter','latex','FontSize',20)



figure (4)
histogram(H_beta_Train(:,1),'normalization','probability','FaceColor','r')
hold on
histogram(H_beta_Test(:,1),'normalization','probability','FaceColor','b')
set(gca,'FontSize',14,'xscale','log')
xlabel('$H(\beta)$','Interpreter','latex','FontSize',24)
ylabel('$f(\cdot)$','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','Interpreter','latex','FontSize',20)
grid on
grid minor

figure (5)
histogram(H_beta_Train(:,33),'normalization','probability','FaceColor','r')
hold on
histogram(H_beta_Test(:,33),'normalization','probability','FaceColor','b')
set(gca,'FontSize',14,'xscale','log')
xlabel('$H(\beta)$','Interpreter','latex','FontSize',24)
ylabel('$f(\cdot)$','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','Interpreter','latex','FontSize',20)
grid on
grid minor

figure (6)
histogram(H_beta_Train(:,60),'normalization','probability','FaceColor','r')
hold on
histogram(H_beta_Test(:,60),'normalization','probability','FaceColor','b')
set(gca,'FontSize',14,'xscale','log')
xlabel('$H(\beta)$','Interpreter','latex','FontSize',24)
ylabel('$f(\cdot)$','Interpreter','latex','FontSize',24)
legend('Train Set','Test Set','Interpreter','latex','FontSize',20)
grid on
grid minor


%% Parameter evolution

clear

load('Results_May24_STUDENT_Test_1_CV_1.mat')


for i = 1:length(VecL)

    beta = M_beta{i};
    beta1(i) = beta(2);
    beta2(i) = beta(3);
    beta3(i) = beta(4);
    beta4(i) = beta(5);
    beta5(i) = beta(6);
    beta6(i) = beta(7);

end



figure
semilogx(VecL,beta1,'linewidth',2)
hold on
semilogx(VecL,beta2,'linewidth',2)
semilogx(VecL,beta3,'linewidth',2)
semilogx(VecL,beta4,'linewidth',2)
semilogx(VecL,beta5,'linewidth',2)
semilogx(VecL,beta6,'linewidth',2)
set(gca,'FontSize',14)
xlabel('$\lambda$','Interpreter','latex','FontSize',24)
ylabel('Values','Interpreter','latex','FontSize',24)
legend('\beta_{1}','\beta_2','\beta_3','\beta_4','\beta_5','\beta_6','Interpreter','tex','FontSize',20)
grid on
grid minor








