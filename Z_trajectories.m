clc
clear


%CODE FOR Z-ENCODED PREFERENCE AND AUTOSOMAL TRAIT
%Z is the focal preference allele, and T is the focal trait allele
%see Muralidhar 2019 (https://www.nature.com/articles/s41586-019-1271-7)
%for more details

% note that the sexually antagonistic trait allele (T) described below is beneficial in males, and
% deleterious in females, in constrast to the trait described in the cases of
% W-linked and X-linked preferences.


% Genotype order: females: zw tt, zw Tt, zw TT, Zw tt, Zw Tt, Zw TT
%                 males:   zz tt, zz Tt, zz TT, Zz tt, Zz Tt, Zz TT, ZZ tt, ZZ Tt, ZZ TT
% all: zz tt, zz Tt, zz TT, Zz tt, Zz Tt, Zz TT, ZZ tt, ZZ Tt, ZZ TT, zw tt, zw Tt, zw TT, Zw tt, Zw Tt, Zw TT

gens = 50000; %number of generations

sf = 0.1; % viability cost of mutant trait allele to homozyg. females 
sm = 0.05; % viability benefit of mutant trait allele to homozyg. males

h = 0.5; % dominance of mutant trait allele, T, with respect to wild-type allele, t
 
alpha = 5; %strength of preference for TT males (see Figure 2 in Muralidhar 2019)

vf = [1; 1-h*sf; 1-sf; 1; 1-h*sf; 1-sf]; % viabilities of female genotypes
vm = [1; 1+h*sm; 1+sm; 1; 1+h*sm; 1+sm; 1; 1+h*sm; 1+sm]; % viabilities of male genotypes

u = 0;      % mutation rate at preference locus, probability per replication
v = 0.001;  % mutation rate at trait locus, probability per replication


%% mate propensity matrix

% used to generate M(i,j), the fraction, out of all matings, that are between an i female and a j male 

Mp = [1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1;
      1,alpha^(h),alpha,1,alpha^(h),alpha,1,alpha^(h),alpha;
      1,alpha^(h),alpha,1,alpha^(h),alpha,1,alpha^(h),alpha;
      1,alpha^(h),alpha,1,alpha^(h),alpha,1,alpha^(h),alpha];
  
%% gamete fractions  

% s(i,k) is the fraction of type i male's gametes that are type k
% sperm types ordered zt, zT, Zt, ZT
  
s = [(1-u)*(1-v), (1-u)*v, u*(1-v), u*v;
       (1-u)/2, (1-u)/2, u/2, u/2;
       (1-u)*v, (1-u)*(1-v), u*v, u*(1-v);
       (1-v)/2, v/2, (1-v)/2, v/2;
       1/4, 1/4, 1/4, 1/4;
       v/2, (1-v)/2, v/2, (1-v)/2;
       u*(1-v), u*v, (1-u)*(1-v), (1-u)*v;
       u/2, u/2, (1-u)/2, (1-u)/2;
       u*v, u*(1-v), (1-u)*v, (1-u)*(1-v)];
   
% e(i,k) is the fraction of type i female's gametes that are type k
% egg types ordered zt, zT, Zt, ZT, wt, wT
   
e = [(1-u)*(1-v)/2, (1-u)*v/2, u*(1-v)/2, u*v/2, (1-v)/2, v/2;
         (1-u)/4, (1-u)/4, u/4, u/4, 1/4, 1/4;
         (1-u)*v/2, (1-u)*(1-v)/2, u*v/2, u*(1-v)/2, v/2, (1-v)/2;
         u*(1-v)/2, u*v/2, (1-u)*(1-v)/2, (1-u)*v/2, (1-v)/2, v/2;
         u/4, u/4, (1-u)/4, (1-u)/4, 1/4, 1/4;
         u*v/2, u*(1-v)/2, (1-u)*v/2, (1-u)*(1-v)/2, v/2, (1-v)/2];
     
%% offspring matrix

% O(i,j,k) is the fraction, out of the offspring produced by (i,j) matings,
% that are type k (among all ordered genotypes)

O = zeros(6,9,15);

for j = 1:6
    for i = 1:9
        
        O(j,i,:) = [s(i,1)*e(j,1), s(i,1)*e(j,2)+s(i,2)*e(j,1), s(i,2)*e(j,2),...
                    s(i,1)*e(j,3)+s(i,3)*e(j,1), s(i,1)*e(j,4)+s(i,2)*e(j,3)+s(i,3)*e(j,2)+s(i,4)*e(j,1), s(i,2)*e(j,4)+s(i,4)*e(j,2),...
                    s(i,3)*e(j,3), s(i,3)*e(j,4)+s(i,4)*e(j,3), s(i,4)*e(j,4),...
                    s(i,1)*e(j,5), s(i,1)*e(j,6)+s(i,2)*e(j,5), s(i,2)*e(j,6),...
                    s(i,3)*e(j,5), s(i,3)*e(j,6)+s(i,4)*e(j,5), s(i,4)*e(j,6)];
                
    end
end

%% loop over generations

p = zeros(6,gens); % female genotype frequencies
q = zeros(9,gens); % male genotype frequencies

startZ = 0.01; % starting frequency of mutant Z allele
startT = 0.01; % starting frequency of mutant T allele


%starting frequencies of all genotypes, in HWE
start_p = [(1-startZ)*(1-startT)^2; (1-startZ)*2*startT*(1-startT); (1-startZ)*startT^2; startZ*(1-startT)^2;...
          startZ*2*startT*(1-startT); startZ*startT^2];
      
start_q = [(1-startZ)^2*(1-startT)^2; (1-startZ)^2*2*startT*(1-startT); (1-startZ)^2*startT^2; ...
          2*startZ*(1-startZ)*(1-startT)^2; 2*startZ*(1-startZ)*2*startT*(1-startT); 2*startZ*(1-startZ)*startT^2; ...
          startZ^2*(1-startT)^2; startZ^2*2*startT*(1-startT); startZ^2*startT^2];
      
%normalize starting frequencies
p(:,1) = start_p/sum(start_p);
q(:,1) = start_q/sum(start_q);


for g = 2:gens
    
    pv = p(:,g-1).*vf; pv = pv/sum(pv); % female genotypes frequencies after viability selection
    qv = q(:,g-1).*vm; qv = qv/sum(qv); % male genotypes frequencies after viability selection
    
    M = repmat(qv',6,1).*Mp;
    M = M./repmat(sum(M')',1,9);
    M = M.*repmat(pv,1,9);
    M = M/sum(sum(M));
    % M(i,j) is the fraction, out of all matings, that are between an i female and a j male
    
    offspring = repmat(M,1,1,15).*O;
    
    geno = sum(squeeze(sum(offspring)));
    
    q(:,g) = geno(1:9)/sum(geno(1:9));
    p(:,g) = geno(10:15)/sum(geno(10:15));
    
end

figure('Position',[100 100 500 400])
set(gcf,'Color','w');
set(gca,'FontName','times');


plot(1:gens, sum(p(4:6,:)), 'Linewidth', 2)
hold on
plot(1:gens, sum(p([2,5],:))/2 + sum(p([3,6],:)), 'Linewidth', 2)
axis([0 gens 0 1])
legend('PreferenceZ', 'Trait')
axis([0 gens 0 1])

xlabel('Generations')
ylabel('Frequency')

