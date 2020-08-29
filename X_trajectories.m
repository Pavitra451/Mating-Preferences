clc
clear


%CODE FOR X-ENCODED PREFERENCE AND AUTOSOMAL TRAIT
%X is the focal preference allele, and T is the focal trait allele
%see Muralidhar 2019 Methods (https://www.nature.com/articles/s41586-019-1271-7) for more details


% Genotype order: females: xx tt, xx Tt, xx TT, Xx tt, Xx Tt, Xx TT, XX tt, XX Tt, XX TT
%                 males:   xy tt, xy Tt, xy TT, Xy tt, Xy Tt, Xy TT
% all: xx tt, xx Tt, xx TT, Xx tt, Xx Tt, Xx TT, XX tt, XX Tt, XX TT, xy tt, xy Tt, xy TT, Xy tt, Xy Tt, Xy TT

gens = 50000; %generations

sf = 0.075; % viability benefit of mutant trait allele to homozyg. females  
sm = 0.1; % viability cost of mutant trait allele to homozyg. males

hT = 0.5; % dominance of mutant trait allele, T, with respect to wild-type allele, t
hX = 0.5; % dominance of mutant preference allele, X, with respect to wild-type allele, x

alpha = 5; % strength of preference of XX females for TT males (see Methods and Figure 2 in Muralidhar 2019)

vf = [1; 1+hT*sf; 1+sf; 1; 1+hT*sf; 1+sf; 1; 1+hT*sf; 1+sf]; % viabilities of female genotypes
vm = [1; 1-hT*sm; 1-sm; 1; 1-hT*sm; 1-sm];                   % viabilities of male genotypes

u = 0;      % mutation rate at preference locus, probability per replication
v = 0.001;  % mutation rate at trait locus, probability per replication


%% mate propensity matrix

% used to generate M(i,j), the fraction, out of all matings, that are between an i female and a j male 

Mp = [1,1,1,1,1,1;
      1,1,1,1,1,1;
      1,1,1,1,1,1;
      1,alpha^(hX*hT),alpha^hX,1,alpha^(hX*hT),alpha^hX;
      1,alpha^(hX*hT),alpha^hX,1,alpha^(hX*hT),alpha^hX;
      1,alpha^(hX*hT),alpha^hX,1,alpha^(hX*hT),alpha^hX;
      1,alpha^hT,alpha,1,alpha^hT,alpha;
      1,alpha^hT,alpha,1,alpha^hT,alpha;
      1,alpha^hT,alpha,1,alpha^hT,alpha];
  
%% gamete fractions  

% e(i,k) is the fraction of type i female's gametes that are type k
% egg types ordered xt, xT, Xt, XT
  
e = [(1-u)*(1-v), (1-u)*v, u*(1-v), u*v;
       (1-u)/2, (1-u)/2, u/2, u/2;
       (1-u)*v, (1-u)*(1-v), u*v, u*(1-v);
       (1-v)/2, v/2, (1-v)/2, v/2;
       1/4, 1/4, 1/4, 1/4;
       v/2, (1-v)/2, v/2, (1-v)/2;
       u*(1-v), u*v, (1-u)*(1-v), (1-u)*v;
       u/2, u/2, (1-u)/2, (1-u)/2;
       u*v, u*(1-v), (1-u)*v, (1-u)*(1-v)];
   
% s(i,k) is the fraction of type i male's gametes that are type k
% sperm types ordered xt, xT, Xt, XT, yt, yT
   
s = [(1-u)*(1-v)/2, (1-u)*v/2, u*(1-v)/2, u*v/2, (1-v)/2, v/2;
         (1-u)/4, (1-u)/4, u/4, u/4, 1/4, 1/4;
         (1-u)*v/2, (1-u)*(1-v)/2, u*v/2, u*(1-v)/2, v/2, (1-v)/2;
         u*(1-v)/2, u*v/2, (1-u)*(1-v)/2, (1-u)*v/2, (1-v)/2, v/2;
         u/4, u/4, (1-u)/4, (1-u)/4, 1/4, 1/4;
         u*v/2, u*(1-v)/2, (1-u)*v/2, (1-u)*(1-v)/2, v/2, (1-v)/2];
     
%% offspring matrix

% O(i,j,k) is the fraction, out of the offspring produced by (i,j) matings, 
% that are type k (among all ordered genotypes)

O = zeros(9,6,15);

for i = 1:9
    for j = 1:6
        
        O(i,j,:) = [e(i,1)*s(j,1), e(i,1)*s(j,2)+e(i,2)*s(j,1), e(i,2)*s(j,2),...
                    e(i,1)*s(j,3)+e(i,3)*s(j,1), e(i,1)*s(j,4)+e(i,2)*s(j,3)+e(i,3)*s(j,2)+e(i,4)*s(j,1), e(i,2)*s(j,4)+e(i,4)*s(j,2),...
                    e(i,3)*s(j,3), e(i,3)*s(j,4)+e(i,4)*s(j,3), e(i,4)*s(j,4),...
                    e(i,1)*s(j,5), e(i,1)*s(j,6)+e(i,2)*s(j,5), e(i,2)*s(j,6),...
                    e(i,3)*s(j,5), e(i,3)*s(j,6)+e(i,4)*s(j,5), e(i,4)*s(j,6)];
                
    end
end

%% loop over generations 

p = zeros(9,gens); % female genotype frequencies
q = zeros(6,gens); % male genotype frequencies

startX = 0.01; % starting frequency of mutant X allele
startT = 0.01; % starting frequency of mutant T allele

%starting frequencies of all genotypes, in HWE
start_p = [(1-startX)^2*(1-startT)^2; (1-startX)^2*2*startT*(1-startT); (1-startX)^2*startT^2; ...
          2*startX*(1-startX)*(1-startT)^2; 2*startX*(1-startX)*2*startT*(1-startT); 2*startX*(1-startX)*startT^2; ...
          startX^2*(1-startT)^2; startX^2*2*startT*(1-startT); startX^2*startT^2];
start_q = [(1-startX)*(1-startT)^2; (1-startX)*2*startT*(1-startT); (1-startX)*startT^2; startX*(1-startT)^2;...
          startX*2*startT*(1-startT); startX*startT^2];
      
%normalize starting frequencies
p(:,1) = start_p/sum(start_p);
q(:,1) = start_q/sum(start_q);


for g = 2:gens
    
    pv = p(:,g-1).*vf; pv = pv/sum(pv); % female genotypes frequencies after viability selection
    qv = q(:,g-1).*vm; qv = qv/sum(qv); % male genotypes frequencies after viability selection
    
    M = repmat(qv',9,1).*Mp;
    M = M./repmat(sum(M')',1,6);
    M = M.*repmat(pv,1,6);
    M = M/sum(sum(M));
    % M(i,j) is the fraction, out of all matings, that are between an i female and a j male
    
    offspring = repmat(M,1,1,15).*O;
    
    geno = sum(squeeze(sum(offspring)));
    
    p(:,g) = geno(1:9)/sum(geno(1:9));
    q(:,g) = geno(10:15)/sum(geno(10:15));
    
end



figure('Position',[100 100 500 400])
set(gcf,'Color','w');
set(gca,'FontName','times');

plot(1:gens, sum(q(4:6,:)), 'Linewidth', 2)
hold on
plot(1:gens, sum(q([2,5],:))/2 + sum(q([3,6],:)), 'Linewidth', 2)
legend('PreferenceX', 'Trait')
axis([0 gens 0 1])

xlabel('Generations')
ylabel('Frequency')

    
    
