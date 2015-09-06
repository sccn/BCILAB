function [ALPHA ALPHA_BAR THETA FIT phi_t]=smooth_mcmc(Y,K,Q,nloop)
     
N=size(Y,1);
T=size(Y{1},1);
time=1:T;

%% set up basis design matrices

knots=quantile(time,0);
order=4;
for(k=1:(K-order+1))
    knots=[knots quantile(time,k/(K-order+1))];
end
knots=knots';

delta=1/1000;
time_cont=(1:(10*T))/10;
T_cont=length(time_cont);
phi_t=spcol(augknt(knots,4),4,time)';
for k=1:K
    if k>1
        for kk=1:(k-1)
            phi_t(k,:)=phi_t(k,:)-((phi_t(k,:)*phi_t(kk,:)')/(phi_t(kk,:)*phi_t(kk,:)'))*phi_t(kk,:);
        end
    end
    phi_t(k,:)=phi_t(k,:)/sqrt(phi_t(k,:)*phi_t(k,:)');
end    
% plot(phi_t(1,:))
% hold on
% for k=2:K
%     plot(phi_t(k,:))
% end
% hold off



    
%% initialize parameters

ALPHA=zeros(Q,N,nloop+1);
ALPHA_BAR=zeros(Q,nloop+1);
THETA=zeros(K,Q,nloop+1);
THETA(:,1,1)=mvnrnd(zeros(K,1),eye(K));
THETA(:,1,1)=THETA(:,1,1)/sqrt(THETA(:,1,1)'*THETA(:,1,1));
for(q=2:Q)
    THETA(:,q,1)=mvnrnd(zeros(K,1),eye(K));
    for(r=1:(q-1))
        THETA(:,q,1)=THETA(:,q,1)-(THETA(:,q,1)'*THETA(:,q-r,1)/...
            (THETA(:,q-r,1)'*THETA(:,q-r,1)))*THETA(:,q-r,1);
    end
    THETA(:,q,1)=THETA(:,q,1)/sqrt(THETA(:,q,1)'*THETA(:,q,1));
end
SIGMA_EPS=.1*ones(1,nloop+1);
FIT=zeros(T,N,nloop+1);


%%Hyperparameters

phi=100;
pi1=0.01;
pi2=0.01;


%% run MCMC to smooth data


for(iter=1:nloop)


    iter


    % Draw ALPHA

    for(i=1:N)
        Sigma_alpha_i=inv(eye(Q)+THETA(:,:,iter)'*phi_t(:,1)*phi_t(:,1)'*THETA(:,:,iter))/SIGMA_EPS(iter);
        mu_alpha_i=ALPHA_BAR(:,iter)+THETA(:,:,iter)'*phi_t(:,1)*Y{i}(1)/SIGMA_EPS(iter);
        for(t=2:T)
            Sigma_alpha_i=Sigma_alpha_i+THETA(:,:,iter)'*phi_t(:,t)*phi_t(:,t)'*THETA(:,:,iter)/SIGMA_EPS(iter);
            mu_alpha_i=mu_alpha_i+THETA(:,:,iter)'*phi_t(:,t)*Y{i}(t)/SIGMA_EPS(iter);
        end
        Sigma_alpha_i=inv(Sigma_alpha_i);
        mu_alpha_i=Sigma_alpha_i*mu_alpha_i;
        ALPHA(:,i,iter+1)=mvnrnd(mu_alpha_i,Sigma_alpha_i)';
    end


    % Draw ALPHA_BAR

    Sigma_alpha_bar=eye(Q)/(N+1/phi);
    mu_alpha_bar=0;
    for(i=1:N)
        mu_alpha_bar=mu_alpha_bar+ALPHA(:,i,iter+1);
    end
    mu_alpha_bar=Sigma_alpha_bar*mu_alpha_bar;
    ALPHA_BAR(:,iter+1)=mvnrnd(mu_alpha_bar,Sigma_alpha_bar)';


    
    % Draw SIGMA_EPS

    pi1_eps=pi1+N*T/2;
    pi2_eps=pi2;
    for(i=1:N)
        pi2_eps=(Y{i}(1)-phi_t(:,1)'*THETA(:,:,iter+1)*ALPHA(:,i,iter+1))^2/2;
        for(t=2:T)
           pi2_eps=pi2_eps+(Y{i}(t)-phi_t(:,t)'*THETA(:,:,iter+1)*ALPHA(:,i,iter+1))^2/2;
        end
    end
    SIGMA_EPS(:,iter+1)=1/gamrnd(pi1_eps,1/pi2_eps);
   

    % Draw THETA

    Sigma_theta=eye(Q*K)/phi;
    mu_theta=zeros(Q*K,1);
    for(i=1:N)
        Alpha_i=ALPHA(:,i,iter+1)';
        for(t=1:T)
            Phi_D_t=kron(eye(Q),phi_t(:,t)');
            Sigma_theta=Sigma_theta+Phi_D_t'*Alpha_i'*Alpha_i*Phi_D_t/SIGMA_EPS(iter+1);
            mu_theta=mu_theta+Phi_D_t'*Alpha_i'*Y{i}(t)/SIGMA_EPS(iter+1);
        end
    end
    Sigma_theta=.5*(Sigma_theta+Sigma_theta');
    Sigma_theta=inv(Sigma_theta);
    Sigma_theta=.5*(Sigma_theta+Sigma_theta');
    mu_theta=Sigma_theta*mu_theta;
    THETA(:,:,iter+1)=reshape(mvnrnd(mu_theta,Sigma_theta)',K,Q);
    Sigma_psi=zeros(K);
    for i=1:N
        Sigma_psi=Sigma_psi+THETA(:,:,iter+1)*(ALPHA(:,i,iter+1)-ALPHA_BAR(:,iter+1))*...
            (ALPHA(:,i,iter+1)-ALPHA_BAR(:,iter+1))'*THETA(:,:,iter+1)'/N;
    end
    %Sigma_psi=THETA(:,:,iter+1)*THETA(:,:,iter+1)';
    %Sigma_psi=.5*(Sigma_psi+Sigma_psi');
    %[V D]=eig(Sigma_psi);
    %THETA(:,:,iter+1)=V(:,(K-Q+1):K)*sqrt(D((K-Q+1):K,(K-Q+1):K));
    %THETA(:,:,iter+1)=V(:,(K-Q+1):K);
            
    
    % Compute Fits
    
    for i=1:N
        FIT(:,i,iter+1)=phi_t'*THETA(:,:,iter+1)*ALPHA(:,i,iter+1);
    end
        
        
        
end


%save smooth_mcmc_outputs 

