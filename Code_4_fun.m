function [] = Code_4_fun(X,y,z,XT,yT,zT,F,VecL,Name)


tStart = tic;

xi_ini=zeros(length(y),1);

n_opt = 50;

Vec_a_ini = logspace(-5,15,n_opt);
Vec_Ephi_ini = logspace(-5,15,n_opt);

[A_ini,Ephi_ini] = meshgrid(Vec_a_ini,Vec_Ephi_ini);

if F==1
    z_fin=z;        
else
    z_fin=zT;
end

Tasa_Exito_Train=zeros(length(VecL),1);
Tasa_Exito_Test=zeros(length(VecL),1);

Dif_Train=zeros(length(VecL),1);
Dif_Test=zeros(length(VecL),1);

p97_5_TE_Train=zeros(length(VecL),1);
p2_5_TE_Train=zeros(length(VecL),1);

p97_5_TE_Test=zeros(length(VecL),1);
p2_5_TE_Test=zeros(length(VecL),1);



n_cambios = -length(VecL);

t_iter=0;

for k=1:n_opt
    tstart=tic;

    for k0=1:n_opt

        a_ini = A_ini(k,k0);
        b_ini = a_ini./Ephi_ini(k,k0); 

        Progress = round([100.*k/n_opt , 100.*k0/n_opt , n_cambios, (n_opt-k)*t_iter/60],0)


    
        for i=1:length(VecL)    
          
        
        l=VecL(i);


        [beta, Sigma, invSigma, logdetV, E_a, L, a_opt, b_ini, xi_ini,vec_a,vec_b] = vb_logit_fit_Sigma_Fin(X,y,z_fin',l,a_ini,b_ini,xi_ini);


        Nsim = 10^3;

        %% Train

        p_y_train = vb_logit_pred(X, beta, Sigma, invSigma);

        M_ypred_train = 2*(repmat(p_y_train,1,Nsim)>=rand(length(y),Nsim))-1;

        M_Successes_train = (M_ypred_train == repmat(y,1,Nsim));

        V_Tasa_Exito_Train = sum(M_Successes_train)./length(y);

        Tasa_Exito_Train_aux = mean(V_Tasa_Exito_Train);
        p97_5_TE_Train_aux = prctile(V_Tasa_Exito_Train,97.5);
        p2_5_TE_Train_aux = prctile(V_Tasa_Exito_Train,2.5);


        %% Test

        p_y_test = vb_logit_pred(XT, beta, Sigma, invSigma);

        M_ypred_test = 2*(repmat(p_y_test,1,Nsim)>=rand(length(yT),Nsim))-1;

        M_Successes_test = (M_ypred_test == repmat(yT,1,Nsim));

        V_Tasa_Exito_Test = sum(M_Successes_test)./length(yT);

        Tasa_Exito_Test_aux = mean(V_Tasa_Exito_Test);
        p97_5_TE_Test_aux = prctile(V_Tasa_Exito_Test,97.5);
        p2_5_TE_Test_aux = prctile(V_Tasa_Exito_Test,2.5);


        if (Tasa_Exito_Train_aux-Tasa_Exito_Train(i)) > 1e-4

            n_cambios = n_cambios +1;

            Tasa_Exito_Train(i) = Tasa_Exito_Train_aux;
            p97_5_TE_Train(i) = p97_5_TE_Train_aux;
            p2_5_TE_Train(i) = p2_5_TE_Train_aux;

            Dif_Train(i) = trace(z'*Sigma*z)+(beta'*z)^2;
            
        
            Tasa_Exito_Test(i) = Tasa_Exito_Test_aux;
            p97_5_TE_Test(i) = p97_5_TE_Test_aux;
            p2_5_TE_Test(i) = p2_5_TE_Test_aux;

            Dif_Test(i) = trace(zT'*Sigma*zT)+(beta'*zT)^2; 


            M_beta{i} = beta;
            M_Sigma{i} = Sigma;
            M_E_a{i} = E_a;
            M_L{i} = L;
            M_vec_a{i} = vec_a';
            M_vec_b{i} = vec_b';            

        end


        end

      
    end

    t_iter = toc(tstart);

end


%%


Lim = length(VecL);

figure 
semilogx(VecL(1:Lim),Tasa_Exito_Train(1:Lim),'r','linewidth',2)
hold on
semilogx(VecL(1:Lim),p97_5_TE_Train(1:Lim),'--r','linewidth',2)
semilogx(VecL(1:Lim),p2_5_TE_Train(1:Lim),'--r','linewidth',2)
semilogx(VecL(1:Lim),Tasa_Exito_Test(1:Lim),'b','linewidth',2)
semilogx(VecL(1:Lim),p97_5_TE_Test(1:Lim),'--b','linewidth',2)
semilogx(VecL(1:Lim),p2_5_TE_Test(1:Lim),'--b','linewidth',2)
grid on
grid minor
set(gca,'FontSize',14)
xlabel('$\lambda$','Interpreter','latex','FontSize',24)
ylabel('Success Rate','Interpreter','latex','FontSize',24)
legend('Train Set','95\% PI (Train)','','Test Set','95\% PI (Test)','Interpreter','latex','FontSize',20)

pause(0.01)     % To plot instantly and not wait for the full loop to finish

tIter = toc(tStart);

save(Name);





end
