function [trial]=DEgenerator(p,objF,minVar,maxVar)
%% input parameters
    % p -- the parent population/the target vectors
    % objF -- the objective function value of the parent population
    % minVar -- the lower bound of the tested fucntion
    % maxVar -- the upper bound of the tested function
%% output parameters
    % trial -- the generated offspring population/the generated trail vectors
    
%%    
% calculate the size of the population p(popsize) and the number of dimensions(n) of
% the tested function
[popsize,n]=size(p);

% assign a matrix to record the mutant vectors 
mutant=zeros(popsize,n);

% assign an array to record the flag which decides the manner of 
% crossover(arithmatic crossover or binomial crossover) 
flag=zeros(popsize,1);

% set the crossover rate for each vector
% the crossover rate is set as 0.1,0.2 or 1.0 with equal probability
% randomly
CR=ones(popsize,1)*0.1;
randArray=rand(popsize,1);
oneIndex= randArray <= 1/3;
CR(oneIndex)=1.0;
point2Index= randArray > 1/3 & randArray <= 2/3;
CR(point2Index)=0.2;


% set the scaling factor for each vector
% the scaling factor is set as 0.6,0.8 or 1.0 with equal probability
% randomly
F=zeros(popsize,1);
randPerm=randperm(popsize)';
F(randPerm(1:round(popsize/3)))=1.0;
F(randPerm(round(popsize/3)+1:round(popsize/3*2)))=0.8;
F(randPerm(round(popsize/3*2)+1:popsize))=0.6;



% the process of crossover for each target vector
for i=1:popsize
    
    % an random value to decide which crossover strategy is to be chosen
    randOper=rand;
    if  randOper < 0.5
       % DE/current-to-best
       % select the best solution which is the solution with minimum
       % objective function value
       [~,minObjFIndex]=min(objF); 
       bestTarget=p(minObjFIndex,:); 
       % randomly select other three distinct vectors which are also different 
       %from both the best vector above and the target vector 
       indexset=1:popsize;
       indexset(i)=[];
       r1=floor(rand*(popsize-1))+1;
       xr1=indexset(r1);
       indexset(r1)=[];
       r2=floor(rand*(popsize-2))+1;
       xr2=indexset(r2);
       indexset(r2)=[];
       r3=floor(rand*(popsize-3))+1;
       xr3=indexset(r3);
       
       % set the ith value of flag as 0 which means that the ith mutant
       % vector would adopt the binomial crossover
       flag(i)=0; 
       
       % utilize the DE/current-to-best strategy to generate the mutant
       % vector
      mutant(i,:)=p(xr1,:)+rand*(bestTarget-p(xr1,:))+F(i)*(p(xr2,:)-p(xr3,:));
    else
         
       % DE/current-to-rand  
       % select other three distinct vectors which are also different from
       % the target vector randomly
       index=[i];
       indexset=1:popsize;
       indexset(index)=[];
       r1=floor(rand*(popsize-1))+1;
       xr1=indexset(r1);
       indexset(r1)=[];
       r2=floor(rand*(popsize-2))+1;
       xr2=indexset(r2);
       indexset(r2)=[];
       r3=floor(rand*(popsize-3))+1;
       xr3=indexset(r3);
       
       % set the ith value of flag as 1 which means that the ith mutant
       % vector would adpot the arithmatic crossover
       flag(i)=1; 
       
       % utilize the DE/current-to-rand strategy to generate the mutant
       % vector
       mutant(i,:)=p(i,:)+rand*(p(xr1,:)-p(i,:))+F(i)*(p(xr2,:)-p(xr3,:));
     end     
end

  %modify the values of some elements violating the boundary constraint,
  %by reflecting them back from the violated boundary
  minMatrix=repmat(minVar,popsize,1); 
  maxMatrix=repmat(maxVar,popsize,1); 

  maxnum=find(mutant>maxMatrix);
  mutant(maxnum)=max((2*maxMatrix(maxnum)-mutant(maxnum)), minMatrix(maxnum));
  minnum=find(mutant< minMatrix);
  mutant(minnum)=min((2* minMatrix(minnum)-mutant(minnum)),maxMatrix(minnum));
  clear num;
  
% crossover operator
 index=find(flag==1);
 term=mutant(index,:);
 t=rand(popsize,n)<repmat(CR,1,n);
 jrand=floor(rand(popsize,1)*n)+1;
 jrand=([1:popsize]'-1)*n+jrand;
 t_T=t';
 t_T(jrand)=1;
 t=t_T';
 t_=1-t;
 trial=mutant.*t+p.*t_;    
 trial(index,:)=term;
end