clc
clear

%CODE FOR AUTOSOMAL PREFERENCE AND AUTOSOMAL TRAIT
%P is the focal preference allele, and T is the focal trait allele
%see Muralidhar 2019 Methods (https://www.nature.com/articles/s41586-019-1271-7)
%for more details

% genotypes pt/pt, Pt/pt, Pt/Pt, pT/pt, PT/pt, pT/Pt, PT/Pt, pT/pT, PT/pT, PT/PT
% gametes pt, Pt, pT, PT
gens = 50000; % number of generations

sf = 0.05; %benefit of trait allele to homozygous females
sm = 0.2; %cost of trait allele to homozygous males


hT = 0.5; % dominance of mutant trait allele, T, with respect to wild-type allele, t
hP = 0.5; % dominance of mutant preference allele, P, with respect to wil-type allele, p

alpha = 5; %strength of preference of PP females for TT males (see Methods and Figure 2 in Muralidhar 2019)

vf = [1; 1; 1; 1+hT*sf; 1+hT*sf; 1+hT*sf; 1+hT*sf; 1+sf; 1+sf; 1+sf]; % viabilities of females of the ten different
% genotypes, with genotypes enumrated in order as above (note the two double heterozygotes PT/pt and pT/Pt)
vm = [1; 1; 1; 1-hT*sm; 1-hT*sm; 1-hT*sm; 1-hT*sm; 1-sm; 1-sm; 1-sm]; % " " males

u = 0;      % mutation rate at preference locus, probability per replication
v = 0.001;  % mutation rate at trait locus, probability per replication

r = 0.5; % recombination rate between trait and preference locus

%% Mate propensity matrix 0 0 0 1 1 1 1 2 2 2

% Mp(i,j) is the "propensity" with which mating pairs form between i females and j males, where genotypes are
% enumerated as above

Mp = [1,1,1,1,1,1,1,1,1,1;
      1,1,1, alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^hP, alpha^hP, alpha^hP;
      1,1,1, alpha^hT, alpha^hT, alpha^hT, alpha^hT, alpha, alpha, alpha;
      1,1,1,1,1,1,1,1,1,1;
      1,1,1, alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^hP, alpha^hP, alpha^hP;
      1,1,1, alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^hP, alpha^hP, alpha^hP;
      1,1,1, alpha^hT, alpha^hT, alpha^hT, alpha^hT, alpha, alpha, alpha;
      1,1,1,1,1,1,1,1,1,1;
      1,1,1, alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^(hP*hT), alpha^hP, alpha^hP, alpha^hP;
      1,1,1, alpha^hT, alpha^hT, alpha^hT, alpha^hT, alpha, alpha, alpha];
  

%% gamete fractions

% g(i,j) is the fraction of genotype i's gametes that are type j, where genotypes and gametes are 
% enumerated as above

g = [(1-u)*(1-v), u*(1-v), (1-u)*v, u*v;
     (1-v)/2, (1-v)/2, v/2, v/2;
     u*(1-v), (1-u)*(1-v), u*v, (1-u)*v;
     (1-u)/2, u/2, (1-u)/2, u/2;
     ((1-r)*(u*v + (1-u)*(1-v)) + r*(u*(1-v) + (1-u)*v))/2, ...
     ((1-r)*((1-u)*v + u*(1-v)) + r*((1-u)*(1-v) + u*v))/2, ...
     ((1-r)*(u*(1-v) + (1-u)*v) + r*(u*v + (1-u)*(1-v)))/2, ...
     ((1-r)*((1-u)*(1-v) + u*v) + r*((1-u)*v + u*(1-v)))/2;
     (r*(u*v + (1-u)*(1-v)) + (1-r)*(u*(1-v) + (1-u)*v))/2, ...
     (r*((1-u)*v + u*(1-v)) + (1-r)*((1-u)*(1-v) + u*v))/2, ...
     (r*(u*(1-v) + (1-u)*v) + (1-r)*(u*v + (1-u)*(1-v)))/2, ...
     (r*((1-u)*(1-v) + u*v) + (1-r)*((1-u)*v + u*(1-v)))/2;
     u/2, (1-u)/2, u/2, (1-u)/2;
     (1-u)*v, u*v, (1-u)*(1-v), u*(1-v);
     v/2, v/2, (1-v)/2, (1-v)/2;
     u*v, (1-u)*v, u*(1-v), (1-u)*(1-v)];
 
%% offspring matrix

O = zeros(10,10,10);

% O(i,j,k) is the fraction of offspring from matings between type-i females and type-j males that are type-k

for i = 1:10
    for j = 1:10
        
        O(i,j,:) = [g(i,1)*g(j,1), g(i,1)*g(j,2) + g(j,1)*g(i,2), g(i,2)*g(j,2), g(i,3)*g(j,1) + g(i,1)*g(j,3), ...
                    g(i,4)*g(j,1) + g(i,1)*g(j,4), g(i,3)*g(j,2) + g(i,2)*g(j,3), g(i,4)*g(j,2) + g(i,2)*g(j,4), ...
                    g(i,3)*g(j,3), g(i,4)*g(j,3) + g(i,3)*g(j,4), g(i,4)*g(j,4)];
                
    end
end

%% loop over generations


startP = 0.01; % starting frequency of mutant P allele
startT = 0.01; % starting frequency of mutant T allele

% starting genotypes frequencies, in HWE
start = [(1-startP)^2*(1-startT)^2; 2*startP*(1-startP)*(1-startT)^2; startP^2*(1-startT)^2; ...
         (1-startP)^2*2*startT*(1-startT); 2*startP*(1-startP)*2*startT*(1-startT)/2; 2*startP*(1-startP)*2*startT*(1-startT)/2; ...
         startP^2*2*startT*(1-startT); (1-startP)^2*startT^2; 2*startP*(1-startP)*startT^2; startP^2*startT^2];
     
p(:,1) = start/sum(start);  % Starting frequencies of female genotypes.
                            
                            % p(i,g) will be the frequency of genotype i in generation g in females, before viability selection
                            
                            

q(:,1) = start/sum(start);  % " " male genotypes

for gen = 2:gens
    
    pv = p(:,gen-1).*vf;   
    pv = pv/sum(pv);       % female genotype frequencies after viavility selection
    
    qv = q(:,gen-1).*vm; qv = qv/sum(qv);   % male genotype frequencies after viavility selection
    
   
    
    M = repmat(qv',10,1).*Mp; 
    M = M./repmat(sum(M')',1,10);         
    M = M.*repmat(pv,1,10);      
    M = M/sum(sum(M));            
    % M(i,j) is the fraction, out of all matings, that are between an i female and a j male                             
    
    offspring = repmat(M,1,1,10).*O;    
                                        
    geno = sum(squeeze(sum(offspring))); 
                                         
    p(:,gen) = geno/sum(geno);  % Next generations female genotype frequencies, normalized to sum to one
    q(:,gen) = geno/sum(geno);  % " "  male

    
end

figure('Position',[100 100 500 400])
set(gcf,'Color','w');
set(gca,'FontName','times');   


plot(1:gens, sum(q([2,5,6,9],:))/2 + sum(q([3,7,10],:)), 'Linewidth', 2)    % plots the frequency of the P allele
                                                                            % across time 
hold on
plot(1:gens, sum(q([4,5,6,7],:))/2 + sum(q([8,9,10],:)), 'Linewidth', 2)
legend('Preference', 'Trait')
axis([0 gens 0 1])


  