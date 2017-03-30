function [ Wij ] = LogNWeightMat( p,mu,sig,N,maxw )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% DEV

p = 0.2;
mu = 1;
sig = 0.5;
N = 1000;
maxw = 15;

%%
Wij = exp(sig.*randn(N)+mu);
%Remove connections >max
bigsyn = Wij>maxw;
while sum(bigsyn(:))>0
    newsyn = exp(sig.*randn(1,sum(bigsyn(:)))+mu);
    Wij(bigsyn) = newsyn;
    bigsyn = Wij>maxw;
end


connects = rand(N)<=p;
connects(diag(diag(true(size(connects)))))=0;

indegree = sum(connects,2)./N;

Wij(~connects) =0;

inweight = sum(Wij,2)./N;

[~,sortinweight] = sort(inweight);
[~,sortindegree] = sort(indegree);

%%
numbins = 25;
weightbins = linspace(-1,log10(maxw),numbins);
weightdist = hist(log10(Wij)',weightbins);

figure
imagesc(weightdist(:,sortindegree))
axis xy


%%
weightplot = [1 1 1; makeColorMap([0 0 0.5],[0.9 0 0])];


figure
subplot(2,2,1)
colormap(weightplot)
imagesc(Wij(sortinweight,sortinweight))
colorbar
xlabel('Neuron j');ylabel('Neuron i')

subplot(4,2,2)
hist(indegree)
xlabel('In Degree (# Inputs)')
ylabel('# Neurons')

subplot(4,2,4)
hist(inweight)
xlabel('Mean Input Connection')
ylabel('# Neurons')

subplot(2,2,4)
plot(indegree,inweight,'.')
xlabel('In Degree');ylabel('Mean Input Connection')

subplot(2,2,3)
hist(Wij(Wij~=0),50)
xlabel('Synaptic Weight (mV)')
ylabel('# Synapses')

end

