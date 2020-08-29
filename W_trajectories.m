clc
clear

%CODE FOR W-ENCODED PREFERENCE AND AUTOSOMAL TRAIT
%W is the focal preference allele, and T is the focal trait allele
%see Muralidhar 2019 Methods(https://www.nature.com/articles/s41586-019-1271-7)
%for more details

% Genotype order females:   w tt, w Tt, w TT, W tt, W Tt, W TT
% Genotype order males:     tt, Tt, TT     
% Genotype order general:   w tt, w Tt, w TT, W tt, W Tt, W TT, tt, Tt, TT

gens = 50000;   % number of generations 

sf = 0.05;  % benefit of mutant trait allele to homozyg. females  
sm = 0.2;   % cost of mutant trait allele to homozyg. males

h = 0.5;      % dominance of mutant trait allele, T, with respect to wild-type allele, t

alpha = 5;    % strength of preference for homozygous males (see Figure 2 of Muralidhar 2019)
beta = alpha^h;    % strength of preference for heterozygous males
c = [1,beta,alpha]; 


vf = [1,1+h*sf,1+sf,1,1+h*sf,1+sf]'; % viabilities of the female genotypes
vm = [1,1-h*sm,1-sm]'; % viabilities of the male genotypes

u = 0.000;  % mutation rate at pref. locus, probability per replication
v = 0.001;  % mutation rate at trait locus, probability per replication

%% set up offspring proportion array

% O(i,j,k) will be the proportion of type k offspring produced in 
% a mating between an i female and a j male

e = [(1/2)*(1-u)*(1-v),(1/2)*(1-u)*v,(1/2)*u*(1-v),(1/2)*u*v,(1/2)*(1-v),(1/2)*v;...
     (1/4)*(1-u),(1/4)*(1-u), (1/4)*u,(1/4)*u,(1/4),(1/4);...
     (1/2)*(1-u)*v,(1/2)*(1-u)*(1-v),(1/2)*u*v,(1/2)*u*(1-v),(1/2)*v,(1/2)*(1-v);...
     (1/2)*u*(1-v), (1/2)*u*v, (1/2)*(1-u)*(1-v), (1/2)*(1-u)*v, (1/2)*(1-v), (1/2)*v;...
     (1/4)*u, (1/4)*u, (1/4)*(1-u), (1/4)*(1-u), (1/4), (1/4);...
     (1/2)*u*v,(1/2)*u*(1-v), (1/2)*(1-u)*v, (1/2)*(1-u)*(1-v), (1/2)*v, (1/2)*(1-v)];  % e(i,k) is the proportion of type k 
                                                                                        % eggs produced by an i female
 
s = [1-v,v;1/2,1/2;v,1-v];  % s(j,k) is the proportion of type k sperm produced by a j male

O = zeros(6,3,9);

for i = 1:6
    for j = 1:3
        
        O(i,j,:) = [e(i,1)*s(j,1), e(i,1)*s(j,2)+e(i,2)*s(j,1), e(i,2)*s(j,2),...
                    e(i,3)*s(j,1), e(i,3)*s(j,2)+e(i,4)*s(j,1), e(i,4)*s(j,2),...
                    e(i,5)*s(j,1), e(i,5)*s(j,2)+e(i,6)*s(j,1), e(i,6)*s(j,2)];
                
    end
end

%% loop over generations

startW = 10^(-4);   % starting frequency of mutant W allele
startT =  7.9 * 10^(-3); % starting frequency of mutant T allele

p = zeros(6,gens); % female genotype frequencies
q = zeros(3,gens); % male genotype frequencies

% Genotype order females:   w tt, w Tt, w TT, W tt, W Tt, W TT
% Genotype order males:     tt, Tt, TT     

%starting frequencies of all genotypes, in HWE
p(:,1) = [(1-startW)*(1-startT)^2;(1-startW)*2*startT*(1-startT);(1-startW)*startT^2;startW*(1-startT)^2;startW*2*startT*(1-startT);startW*startT^2];
q(:,1) = [((1-startT)^2);2*startT*(1-startT);startT^2];


for t = 2:gens
    
    pv = (p(:,t-1).*vf)/sum(p(:,t-1).*vf);  % female frequencies post viability selection
    qv = (q(:,t-1).*vm)/sum(q(:,t-1).*vm);  % male frequencies post viability selection
    
    a = c*qv;   % a weighting for the matings of preferenced females    
    pref = [ones(3,3);c/a;c/a;c/a];
    
    M = (pv*qv').*pref; % M(i,j) is the fraction, out of all matings, that are between an i female and a j male
    
    x = sum(sum(repmat(M,1,1,9).*O));
    
    p(:,t) = 2*x(1:6); q(:,t) = 2*x(7:9);
    p(:,t) = p(:,t)/sum(p(:,t)); q(:,t) = q(:,t)/sum(q(:,t));
    
end


figure('Position',[100 100 500 400])
set(gcf,'Color','w');
set(gca,'FontName','times');

plot(1:gens,sum(p(4:6,:)), 'LineWidth', 1)
hold on
plot(1:gens,q(2,:)/2 + q(3,:), 'LineWidth', 1)
legend('PreferenceW', 'Trait')
axis([0 gens 0 1])

xlabel('Generations')
ylabel('Frequency')



