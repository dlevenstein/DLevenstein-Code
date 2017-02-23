

%close all;
clear;clc;


clear;clc;close all
T = 30; % 30 seconds
delta = 0.001; % 10 ms  % real system --> 1 ms
K = T/delta;


% simulate the hidden state S={-1,+1}
S = zeros(1,K);
if rand(1)>0.5
    S(1) = 1;
else
    S(1) = -1;
end

S(1) = 1;
k = 1;
while k<K+2
    if S(k)==-1
        a=log(0.67);b=sqrt(2*(log(0.96)-a));
        u = lognrnd(a,b);
        if u<0.15 | u>1
            u = 0.96;
        end
        L = ceil(u/delta);
        S(k+1:k+L) = ones(1,L);
        k = k+L;
    elseif S(k)==1
        a=log(0.14);b=sqrt(2*(log(0.17)-a));
        u = lognrnd(a,b);
        if u<0.05 | u>1
            u = 0.17;
        end
        L = ceil(u/delta);
        S(k+1:k+L) = -ones(1,L);
        k = k+L;
    end
end
S = S(2:K+1);
figure(2);subplot(411);
stairs([1:K]*delta, 0.5+S/2,'k','linewidth',2);xlim([0.1 8]);%ylim([-1.5 1.5]);
ylabel('state S(t)','fontsize',18);
title('(a)','fontsize',18);



% simulating CIF
% lambda = exp(alpha * count + beta * S)
% take 100 ms history (10 * 10)
C = 4;  % no. of tectrodes


n = zeros(C,K);
lambda = zeros(C,K);
N = zeros(C,K);
for k=1:100
    if S(k) == 1
        %n(:,k) = unidrnd(4,[C 1]) + 3;  % tectrode spike count (not per cell)
        %n(:,k) = ones(C,1);
        tem = rand(C,1);
        ind = find(tem>0.7);
        N(ind,k) = 1;
    elseif S(k) == -1
        tem = rand(C,1);
        ind = find(tem<0.02);
        if length(ind) > 0
            N(ind,k) = 1;
        end
    end
end


 
 alpha = [0.06; 0.05; 0.03; 0.05*ones(C-3,1)]; beta = [3.5; 4; 3.8*ones(C-2,1)];
 

 
for k=101:K
    %lambda(:,k) = exp( alpha .* sum(n(:,k-100:k-1),2) + beta * S(k));
    lambda(:,k) = exp( alpha/4 * sum(sum(N(:,k-100:k-1))) + beta * S(k));


    % simulating Poisson spike trains (for C tectrodes)
    p = lambda(:,k) * delta;
    u = rand(C,1);
    ind = find(u<p);
    N(ind,k) = 1;

end
count = sum(N,1); % total counts across tectrodes



figure(2);
subplot(412);
plot([101:K]*delta,lambda(1,101:K),'r');
hold on;
plot([101:K]*delta,lambda(2,101:K),'b'); 
 plot([101:K]*delta,lambda(3,101:K),'g');
plot([101:K]*delta,lambda(4,101:K),'m');
 
xlim([0.1 8]);ylabel('CIF \lambda(t)','fontsize',18);
title('(b)','fontsize',18);




trains = cell(C,1);
figure(2);
subplot(413);
for c=1:C
    ind = find(N(c,:));
     
  %plot raster
    for i=1:length(ind)
        hold on;plot([ind(i), ind(i)]*delta, [c-.4 c+.4], 'k', 'linewidth',0.5);
    end
end
xlim([0.1 8]);ylim([0 C+1]);
ylabel('MUA','fontsize',18)
title('(c)','fontsize',18);

figure(2);
subplot(414);
for i=1:(K/10)
    Y(i) = sum(count((i-1)*10+1:i*10));  % 10 ms window spike count 
end
rate = Y/0.01/C;
Time = [1:10:K]*delta;
fill([Time; flipdim(Time,1)], [zeros(size(Time)); rate], [3 3 3]/8, 'EdgeColor', [3 3 3]/8);




% smoothing
sigma = 0.5; % 6xsigma --> 30 ms
a = [-1*sigma: sigma: 1*sigma];
kernel =  1/sqrt(2*pi*sigma^2) * exp(-0.5*(a/sigma).^2);
kernel = kernel/sum(kernel);

rate = conv2(rate, kernel, 'same');
hold on;plot([1:10:K]*delta,rate,'g','linewidth',2);

xlim([0.1 8]);ylim([0 200]);xlabel('Time (s)','fontsize',18);
ylabel('Spikes / s','fontsize',18);
title('(d)','fontsize',18);


%%----------------------- Decoding --------------- 

% examination of firing rate across cells per 10 ms ---> HMM


trueS = S(1:10:K);
ncount = zeros(1,K/10);
for t=11:K/10
    ncount(t) = sum(Y(t-10:t-1));  % 100 ms  % spike history
end




figure(1);
n0 = 11; n = K/10; delta = 0.01;
subplot(311);stairs([n0:n]*delta,Y(n0:n)); xlim([0.1 8]);
subplot(312);stairs([n0:n]*delta,trueS(n0:n)); xlim([0.1 8]);ylim([-1.1 1.1])




% initialization

Epoch = 100;
%A = rand(2,2); A0 = A./repmat(sum(A,2),1,2);
A0 = [0.1, 0.9; 0.05, 0.95];

n0 = 11; n = K/10;
[para] = glmfit([0.5+0.5*trueS(n0:n); ncount(n0:n)]', Y(n0:n)', 'poisson', 'link', 'log', 'constant', 'on');
init_mu = para(1)
init_alpha = para(2)
init_beta = para(3)

disp('press any key to pursue decoding')
%pause

[EstS, ProbS, Est_alpha, Est_beta, Est_mu, A, B, p0, loglikeli, LogP] = decode_Poisson(Y(n0:n),ncount(n0:n),A0, init_alpha,init_beta,init_mu,Epoch);

subplot(313);
stairs([n0:n]*delta,EstS,'k-'); xlim([0.1 8]);ylim([-0.05 1.05])
hold on;
%stairs([n0:n]*delta,ProbS(2,:),'k:','linewidth',1); 
%area([n0:n]*delta, ProbS(2,:));
xlim([0.1 8]);ylim([-0.05 1.05])
ylabel('Prob(S)','fontsize',18)
xlabel('Time (s)','fontsize',18)


errorrate = length(find((2*EstS-1)-trueS(n0:n)))/ (n-n0+1)






