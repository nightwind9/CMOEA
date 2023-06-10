function  [p,objF,conV]=Debselect(p,objF,conV,trial,objFtrial,conVtrial)
%%
% feasible solutions are better than infeasible solutions and infeasible
  % solution with less degree of constraint violation is better than the one 
  % with larger degree of constraint violation 
   betterIndex=find(conVtrial < conV);
   p(betterIndex,:)=trial(betterIndex,:);
   objF(betterIndex)=objFtrial(betterIndex);
   conV(betterIndex)=conVtrial(betterIndex);
   
  % feasible solution with less objective function value is better than the
  % one with larger objective function value
   betterIndex=find(conVtrial==conV & objFtrial < objF);      
   p(betterIndex,:)=trial(betterIndex,:);
   objF(betterIndex)=objFtrial(betterIndex);
   conV(betterIndex)=conVtrial(betterIndex);
   
   



