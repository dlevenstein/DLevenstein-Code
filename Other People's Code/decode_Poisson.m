function [S, ProbS, alphaa, betaa, mu, A, B, p0, loglikeli,LogP] = decode_Poisson(Y,ncount,A,alphaa,betaa,mu,Epoch,opt)
%This is code to identify UP/DOWN states using a HMM, from Zhe (sage) Chen 
%based on Chen et al 2009, Neural Computation.
%
%INPUTS
%   S: two-state Markov chain (0/1) or (-1/1)
%   Y: time series of discrete counts 
%   A: transition prob=
%%


if nargin<8; opt = 1; end  % display results
if nargin<7; Epoch = 100; end
if nargin<6; mu = 0; end




d = length(betaa);


n = length(Y);
State = [-1, 1]; 
State = [0, 1];

lambda = zeros(2,n);
for k=1:n
    lambda(1,k) = exp(mu + alphaa * State(1) + betaa' * ncount(:,k));
    lambda(2,k) = exp(mu + alphaa * State(2) + betaa' * ncount(:,k));
end

for i=1:2
    B(i,:) = exp(-lambda(i,:)).*(lambda(i,:).^Y)./factorial(Y);
end


p0 = [0.5; 0.5];




Alpha = zeros(2,n);
Beta = zeros(2,n);
Gamma = zeros(2,n);
Zeta = zeros(2,2,n);


t = 1; Dif_LogP = 10;
while t<Epoch & Dif_LogP > 0;
    % E-step: forward algorithm
    % compute Alpha
    C = zeros(1,n);     % scaling vector to avoid numerical inaccuracies.
    Alpha(:,1) = p0 .* B(:,1);
    C(1) = sum(Alpha(:,1));     
    Alpha(:,1) = Alpha(:,1)/C(1); % scaling    
    for k=1:n-1
        % compute Alpha
        Alpha(:,k+1) = ((Alpha(:,k))'*A)' .* B(:,k+1);
        C(k+1) = sum(Alpha(:,k+1));
        if C(k+1) > 0
            Alpha(:,k+1) = Alpha(:,k+1)/C(k+1);
        else            
            [t,k]          
            error('numerical problem');
        end
    end
    %min(C(n0:n))
    %pause;
    LogP(t) = sum(log(C(1:n)+eps));  % log(0) = NaN


    % E-step: backward algorithm
    % compute Beta
    Beta(:,n) = ones(2,1);
    Beta(:,n) = Beta(:,n)./C(n);
    for k=n-1:-1:1
        % compute Beta
        Beta(:,k) = A * (Beta(:,k+1).*B(:,k+1));
        Beta(:,k) = Beta(:,k)/C(k);
    end
    
    for k=1:n-1
        temp = Alpha(:,k) * (Beta(:,k+1).*B(:,k+1))' .* A;
        if(sum(temp(:))>0) % avoid divide by zero
            Zeta(:,:,k) = temp ./ (sum(temp(:)));
        else
            error('numerical problem');
        end
        Gamma(:,k) = sum(Zeta(:,:,k),2);
    end  
    
    %% compute sufficient statistics %%  
    Z = State(1)*Gamma(1,:) + State(2)*Gamma(2,:);   % Z = E[S]
    %ZZ = E[SS]
    %figure;plot(Z(n0:n));pause;
    
    
    
    % M-step:
    
    % update transition matrix
    p0 = Gamma(:,1);
    temp1 = sum(Zeta,3);
    temp2 = sum(Gamma,2);
    A = temp1./repmat(temp2,1,2);
    
    
    % update alpha and beta      
    option = 1;
    if option == 1       % iterative method (no closed-form solution)
        Iter = 30;
        [alphaa,betaa,mu] = mynewtoncode(Y, Z,ncount,alphaa, betaa, mu, Iter);
    elseif option == 2
        X = [Z; ncount]; 
        [para, dev, stats] = glmfit(X', Y', 'poisson', 'link', 'log', 'constant', 'off');
        alphaa = para(1);
        betaa = para(2);
    elseif option == 3
        X = [Z; ncount]; 
        [para, dev, stats] = glmfit(X', Y', 'poisson', 'link', 'log', 'constant', 'on');
        mu = para(1);
        alphaa = para(2);
        betaa = para(3);    
    else % for debugging
        alphaa = alphaa; betaa = betaa;
    end
    
    
    
    
    for k=1:n
        lambda(1,k) = exp(mu + alphaa * State(1) + betaa' * ncount(:,k));  
        lambda(2,k) = exp(mu + alphaa * State(2) + betaa' * ncount(:,k));  
    end
    for i=1:2
        B(i,:) = exp(-lambda(i,:)).*(lambda(i,:).^Y)./factorial(Y);
    end
    
    
    if opt == 1
        disp(['Epoch ', num2str(t), ', log-likelihood= ', num2str(LogP(t)), ', mu=', num2str(mu), ', alpha=', num2str(alphaa)]); % , ', beta=', num2str(betaa)]);
    end
    
    if t>1;
        Dif_LogP = LogP(t) - LogP(t-1);
    end
    t = t+1;

end  % end EM loop


ProbS = Gamma; 
meanS = State(1)*Gamma(1,:) + State(2)*Gamma(2,:);  % expected mean 
p0 = Gamma(:,1);



% Viterbi algorithm for decoding most likely states
delta = zeros(2,n);
psi = zeros(2,n);
S = zeros(1,n);


% working in log-domain 
p0 = log(p0+eps);
AA = log(A+eps);
BB = log(B+eps);

delta(:,1) = p0 + BB(:,1);   % 2-by-1
psi(:,1) = 0;

for t = 2:n
    [delta(:,t), psi(:,t)] = max( (delta(:,t-1) * ones(1,2) + AA)', [], 2); % maximum of every row
    delta(:,t) = delta(:,t) + BB(:,t);
end



[loglikeli,S(n)] = max(delta(:,n));   % max(2-by-1) 
loglikeli = loglikeli / n;

for t = n-1:-1:1
    S(t) = psi(S(t+1),t+1);
end

p0 = exp(p0);
if min(State) == -1
    S = 2*S-3; % change S{1,2} to S{-1,1}
elseif min(State) == 0
    S = S-1;   % change S{1,2} to S{0,1}  
end











% ------------- subrountine ------------- %


function [alpha_new,beta_new, mu_new] = mynewtoncode(Y, Z, count, alpha_old, beta_old, mu_old, Iter);

mu_new = 0;
alpha_new = 0;
beta_new = zeros(size(beta_old)); 

if size(count,1) > size(count,2); 
    count = count'; 
end
if size(Y,1) > size(Y,2); 
    Y = Y'; 
end



tem0 = sum(Y);
update = mu_old + 0.01*randn(1);
for i = 1:Iter
    g = sum(exp(update + alpha_old * Z + beta_old' * count)) - tem0;
    gprime = sum( exp(update + alpha_old * Z + beta_old' * count)); % derivative w.r.t. update
    update = update - g/gprime;

    mu_new = update;
end



tem1 = sum(Y .* Z);
update = alpha_old + 0.01*randn(1); 
for i = 1:Iter
    g = sum(Z .* exp(mu_new + update * Z + beta_old' * count)) - tem1; 
    gprime = sum(Z.^2 .* exp(mu_new + update * Z + beta_old' * count)); % derivative w.r.t. update 
    update = update - g/gprime;  
    
    alpha_new = update;
end




d = size(count,1); 
if d == 1
    tem2 = sum(Y .* count);
    update = beta_old + 0.01*randn(1);
    for i = 1:Iter
        g = sum(count .* exp(mu_new + alpha_new * Z + update * count)) - tem2;
        gprime = sum(count.^2 .* exp(mu_new + alpha_new * Z + update * count)); % derivative w.r.t. update
        update = update - g/gprime;

        beta_new = update;
    end
elseif d > 1
    %[d,ydim] = size(count);
    %[1,ydim] = size(Y);
    tem2 = sum(repmat(Y,d,1) .* count, 2);
    update = beta_old + 0.01*randn(d,1);
    for i = 1:Iter
        g = sum(count .* repmat( exp(mu_new + alpha_new * Z + update' * count), d, 1), 2) - tem2;  % vector
        gprime = count .* repmat( exp(mu_new + alpha_new * Z + update' * count), d, 1) * count';   % matrix
        update = update - inv(gprime) * g;

        beta_new = update;
    end

end



% tem3 = sum(Y);
% update = mu_old + 0.01*randn(1);
% for i = 1:Iter
%     g = sum(exp(update + alpha_new * Z + beta_new' * count)) - tem3;
%     gprime = sum( exp(update + alpha_new * Z + beta_new' * count)); % derivative w.r.t. update
%     update = update - g/gprime;
% 
%     mu_new = update;
% end

    
    





