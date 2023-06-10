%Please cite the paper: Z.-Z. Liu and Y. Wang. Handling constrained multiobjective %optimization problems with constraints in both the decision and objective spaces. IEEE Transactions on %Evolutionary Computation, in press. DOI: 10.1109/TEVC.2019.2894743  

% Remark: for the reference point of DOC9 is set to [0.1 1.1 1.1].

% If you have any question, please contact us via ywang@csu.edu.cn and zhizhongliu@csu.edu.cn 
clc;clear;tic;
format long;format compact;

problem_set    = {'DOC1','DOC2','DOC3','DOC4','DOC5','DOC6','DOC7','DOC8','DOC9'}; % The number of test problems in this type
decision_dim   = [6 16 10 8 8 11 11 10 11];
itrNum         = 1; % This program is run for 100 times

for p_id = 1:9 % 1:length(problem_set)
    
    itrCounter = 1;
    
    while itrCounter<=itrNum
        if  p_id <=7
          max_gen    = 2000;    % maximum generation number
        else
          max_gen = floor(400000/300) + 1;
        end  
        if p_id <= 7
           popsize    = 100;    % population size
        else
            popsize = 300;
        end
        n          = decision_dim(p_id);     % the dimenstion of the decision vector
        
        pcrossover = 0.9;    % crossover probability
        pmutation  = 1/n;    % mutation probability
        di         = 20;     % crossover distribution index
        dim        = 20;     % mutation distribution index
        gen        = 1;      % the first generation
        lu         = decision_range(problem_set{p_id},n)';
        fit_cv     = cmop_test(problem_set{p_id});
        
     %%  %%%%%%%%%%%%% the first pahse: single objective method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = ones(popsize, 1) * lu(1, :) + rand(popsize, n) .* (ones(popsize, 1) * (lu(2, :) - lu(1, :)));
        [fit,cv] = fit_cv(x');
        m  = size(fit,1); % the nubmer of objectives
        % calculate the overall constraint violation
        cv(cv > 0) = 0;cv = abs(cv);
        cv = sum(cv,1);cv = cv';
        fit_a = fit;
        fit= sum (fit,1); fit = fit';
        
        pf = 0; %the feasibility proportion 
        delta = 1; % the normalized difference
        sita = 1/3;
         
   while  (delta > 0.2 ||pf <sita) && gen<=max_gen*0.9 % the last condition to make sure the algorithm can enter the second phase  
  
            trial=DEgenerator(x,fit,lu(1,:),lu(2,:));
            [fittrial,cvtrial] = fit_cv(trial');
         
            cvtrial(cvtrial > 0) = 0;cvtrial= abs(cvtrial);
            cvtrial = sum(cvtrial,1);cvtrial = cvtrial';
            fittrial= sum (fittrial,1); fittrial = fittrial';
              
            x_all = [x;trial];
              
            [x,fit,cv] = Debselect(x,fit,cv,trial,fittrial,cvtrial); % feasibility rule
            gen = gen + 1;
             
          %% compute the feasiblity proportion pf      
            index_1 = find(cv==0); % find feasible solutions 
            pf = length(index_1)/popsize;
              
          %% nomalize the objectives
       if length(index_1)> 0
              x_new=x(index_1,:);
              [fit_a,cv_a] = fit_cv(x_new');
              fit_max = max(fit_a');
              fit_min = min(fit_a');
              control = 1;
              if control==1
                 fit_max_history = fit_max;
                 fit_min_history = fit_min;
                 control = 0;
              else
                 fit_max_history = max([fit_max_history;fit_max]);
                 fit_min_history = min([fit_min_history;fit_min]);
              end
              
             fit_norm = [];
             fit_a = fit_a';
             for i = 1: size(fit_a,1)
                 fit_norm(i,:) = (fit_a(i,:) - fit_min_history)./(fit_max_history-fit_min_history);    
             end
             fit_norm_sum = sum(fit_norm,2);
        end
           
         %% compute the normalized difference between solutions
               PL =floor(popsize*sita)+1; % pedifined length
             if length(index_1) > PL
               [~,index_2] = sort(fit_norm_sum);
                delta = fit_norm_sum(index_2(PL)) - fit_norm_sum(index_2(1));
             end
         
   end
       
        
   %% %%%%%%%%% the second phase: NSGA-II %%%%%%%%%%%%%%%%%%%%%%%%
    
        % initialize the population
        x = x_all;
       % x = ones(2 * popsize, 1) * lu(1, :) + rand(2 * popsize, n) .* (ones(2 * popsize, 1) * (lu(2, :) - lu(1, :)));
        % evaluate the population
        [fit,cv] = fit_cv(x');
        m        = size(fit,1); % the nubmer of objectives
        % calculate the overall constraint violation
        cv(cv > 0) = 0;cv = abs(cv);
       % temp = max(cv,[],2);temp = max(temp,1e-6*ones(size(temp)));
       % cv = cv ./ repmat(temp,[1,size(cv,2)]);
        cv = sum(cv,1);cv = cv';fit = fit';
        fit_obj =[fit,cv];
        
        [x,fit,fit_obj] = TSR_sort(x,fit,fit_obj,cv);
        gen = gen + 1;
        
        
        % evolutionary process
        while gen<=max_gen
            if mod(gen,50) ==0
                gen
            end
            
            pool_size = popsize;
            % tournament selection operator
            mate_pool_chromosome = tournament_selection(x, fit, pool_size, n, m+1, popsize);
            
            % genetic operator(simulated binary crossover && mutation)
            offspring_chromosome = genetic_operator(mate_pool_chromosome, pool_size, n, lu, pcrossover, pmutation, di, dim);
            
            % Mix the parent population and offspring population,size(x)=popsize*n,size(fit)=popsize*( m+1 )
            x(popsize+1 : 2*popsize, :) = offspring_chromosome;
            
            % Calculate the fitness value of the mixed population
            [fit,cv] = fit_cv(x');
            m        = size(fit,1); % the nubmer of objectives
            % calculate the overall constraint violation
            cv(cv > 0) = 0;cv = abs(cv);
           % temp = max(cv,[],2);temp = max(temp,1e-6*ones(size(temp)));
          %  cv = cv ./ repmat(temp,[1,size(cv,2)]);
            cv = sum(cv,1);cv = cv';fit = fit';
            fit_obj =[fit,cv];
            % Sort and select
            [x,fit,fit_obj] = TSR_sort(x,fit,fit_obj,cv);    
            gen = gen + 1;
        end
        
        % save the final results
        file_path = strcat(pwd,'/TOP_NSGA2_CDP_Results/');
        if ~isdir(file_path)
            mkdir(file_path);
        end
        
        % show the percentage of the running
        fprintf('running percentage = %.3f\n',(itrNum * (p_id - 1.0) + itrCounter) / (itrNum * length(problem_set)));
        
        
        % objective and constraints fit_obj
        objName = strcat(file_path,problem_set{p_id},'_obj_cv_',num2str(itrCounter),'.dat');
        save(objName,'fit_obj','-ascii');
        decName = strcat(file_path,problem_set{p_id},'_dec_',num2str(itrCounter),'.dat');
        save(decName,'x','-ascii');
        itrCounter =itrCounter+1;
    end
    
       index = find(fit_obj(:,end)==0);
       fit_obj=fit_obj(index,:);
       if size(fit_obj,2)==3
        scatter(getcolumn(fit_obj(:,1:2),1),getcolumn(fit_obj(:,1:2),2));
       else
        plot3(fit_obj(:,1),fit_obj(:,2),fit_obj(:,3),'o');
       end
      
 end
toc;