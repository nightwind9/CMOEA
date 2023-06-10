
function range = decision_range(name,dim)

range = ones(dim,2);

switch name
    case 'DOC1'
        range = [0 78 33 27 27 27; 1 102 45 45 45 45]';  
    case 'DOC2'
        range =[0 zeros(1, 15);  1 10 * ones(1, 15)]'; 
    case 'DOC3'
        range =[0 0 0 0 0 0 0 0 0 0.01;  1 1 300 100 200 100 1 100 200 0.03]';
    case 'DOC4'
       range =  [0 -10 -10 -10 -10 -10 -10 -10; 1 10 10 10 10 10 10 10]';  
    case 'DOC5'
        range = [0 0 0 0 100 6.3 5.9 4.5; 1 1000 40 40 300 6.7 6.4 6.25]'; 
    case 'DOC6'
        range = [0 -10 * ones(1, 10); 1 10 * ones(1, 10)]'; 
    case 'DOC7'
        range = [0 zeros(1, 10); 1 10 * ones(1, 10)]'; 
    case 'DOC8'
       range = [0 0 500 1000 5000 100 100 100 100 100; 1 1 1000 2000 6000 500 500 500 500 500]';
    case 'DOC9'
       range = [0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1; 1 1 10 10 10 10 10 10 10 10 10]';      
end
end