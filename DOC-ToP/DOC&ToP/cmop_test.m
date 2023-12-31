
function fobj = cmop_test(name)

switch name
    case  'DOC1'
        fobj= @L1;
    case  'DOC2'
        fobj= @L2;
    case  'DOC3'
        fobj= @L3;
    case  'DOC4'
        fobj= @L4; 
    case  'DOC5'
        fobj= @L5; 
    case  'DOC6'
        fobj= @L6; 
    case 'DOC7'
        fobj= @L7; 
    case 'DOC8'
        fobj= @L8; 
    case 'DOC9'
        fobj= @L9; 
    otherwise
        error('The optimized problem is not exsit!');
end
end




function [y,c] = L1(x)
    x = x';

    % basic multi-objective problem
    g_temp = 5.3578547 * x(:, 4).^2 + 0.8356891 * x(:, 2).* x(:, 6) + 37.293239 * x(:, 2) - 40792.141;
    g = g_temp+30665.5386717834 +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1-sqrt(y(:,1))./g);
  
   %constraints in objective space
    c(:,1) = max( -(y(:,1).^2 + y(:,2).^2-1), 0);
    
   %constraints in decision space
    c(:, 2) = + 85.334407 + 0.0056858 * x(:, 3).* x(:, 6) + 0.0006262 * x(:, 2).* x(:, 5) - 0.0022053 * x(:, 4).* x(:, 6) - 92;
    c(:, 3) = -85.334407 - 0.0056858 * x(:, 3).* x(:, 6) - 0.0006262 * x(:, 2).* x(:, 5) + 0.0022053 * x(:, 4).* x(:, 6);
    c(:, 4) = + 80.51249 + 0.0071317 * x(:, 3).* x(:, 6) + 0.0029955 * x(:, 2).* x(:, 3) + 0.0021813 * x(:, 4).^2 - 110;
    c(:, 5) = -80.51249 - 0.0071317 * x(:, 3).* x(:, 6) - 0.0029955 * x(:, 2).* x(:, 3) - 0.0021813 * x(:, 4).^2 + 90;
    c(:, 6) = + 9.300961 + 0.0047026 * x(:, 4).* x(:, 6) + 0.0012547 * x(:, 2).* x(:, 4) + 0.0019085 * x(:, 4) .* x(:, 5) - 25;
    c(:, 7) = -9.300961 - 0.0047026 * x(:, 4).* x(:, 6) - 0.0012547 * x(:, 2).* x(:, 4) - 0.0019085 * x(:, 4) .* x(:, 5) + 20;
  
    y = y';
    c= -c';
end


function [y,c] = L2(x)

        popsize = size(x,2);
        x = x';
       
        a = [-16 2 0 1 0;
            0 -2 0 0.4 2;
            -3.5 0 2 0 0;
            0 -2 0 -4 -1;
            0 -9 -2 1 -2.8;
            2 0 -4 0 0;
            -1 -1 -1 -1 -1;
            -1 -2 -3 -2 -1;
            1 2 3 4 5;
            1 1 1 1 1];
        b = [-40 -2 -0.25 -4 -4 -1 -40 -60 5 1];
        c1 = [30 -20 -10 32 -10;
            -20 39 -6 -31 32;
            -10 -6 10 -6 -10;
            32 -31 -6 39 -20;
            -10 32 -10 -20 30];
        d = [4 8 10 6 2];
        e = [-15 -27 -36 -18 -12];
        
    % basic multi-objective problem
    g_temp = sum(repmat(c1(1:5, 1)', popsize, 1).* x(:, 12:16), 2).* x(:, 12) + sum(repmat(c1(1:5, 2)', popsize, 1).* x(:, 12:16), 2).* x(:, 13)...
            + sum(repmat(c1(1:5, 3)', popsize, 1).* x(:, 12:16), 2).* x(:, 14) + sum(repmat(c1(1:5, 4)', popsize, 1).* x(:, 12:16), 2).* x(:, 15)...
            + sum(repmat(c1(1:5, 5)', popsize, 1).* x(:, 12:16), 2).* x(:, 16) + 2 * sum(repmat(d, popsize, 1).* x(:, 12:16).^3, 2)...
            - sum(repmat(b, popsize, 1).* x(:, 2:11), 2);
    g = (g_temp-32.6555929502) +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1- (y(:,1)).^(1/3)./g);  
    
    %constraints in objective space
    c(:,1) =max( -(sqrt(y(:,1)) + y(:,2)-1), 0);
        d1(:,1) = max(((y(:,1)-1/8).^2 + (y(:,2) -1+sqrt(1/8) ).^2 - 0.15*0.15),0);
        d1(:,2) = max( ((y(:,1)-1/2).^2 + (y(:,2) -1+sqrt(1/2)).^2 - 0.15*0.15),0);  
        d1(:,3) = max(((y(:,1)-7/8).^2 + (y(:,2) - 1+sqrt(7/8)).^2 - 0.15*0.15),0);    
    c(:,2) = min(d1, [], 2);
    
    % constraints in decision space
    c(:, 3) = -2 * sum(repmat(c1(1:5, 1)', popsize, 1).* x(:, 12:16), 2) - 3 * d(1).* x(:, 12).^2 - e(1) + sum(repmat(a(1:10, 1)', popsize, 1).* x(:, 2:11), 2);
    c(:, 4) = -2 * sum(repmat(c1(1:5, 2)', popsize, 1).* x(:, 12:16), 2) - 3 * d(2).* x(:, 13).^2 - e(2) + sum(repmat(a(1:10, 2)', popsize, 1).* x(:, 2:11), 2);
    c(:, 5) = -2 * sum(repmat(c1(1:5, 3)', popsize, 1).* x(:, 12:16), 2) - 3 * d(3).* x(:, 14).^2 - e(3) + sum(repmat(a(1:10, 3)', popsize, 1).* x(:, 2:11), 2);
    c(:, 6) = -2 * sum(repmat(c1(1:5, 4)', popsize, 1).* x(:, 12:16), 2) - 3 * d(4).* x(:, 15).^2 - e(4) + sum(repmat(a(1:10, 4)', popsize, 1).* x(:, 2:11), 2);
    c(:, 7) = -2 * sum(repmat(c1(1:5, 5)', popsize, 1).* x(:, 12:16), 2) - 3 * d(5).* x(:, 16).^2 - e(5) + sum(repmat(a(1:10, 5)', popsize, 1).* x(:, 2:11), 2);
    y = y';
    c= -c';
end

function [y,c] = L3(x)

     x = x';
    % basic multi-objective problem
    g_temp = -9.* x(:, 6) - 15.* x(:, 9) + 6.* x(:, 2) + 16.* x(:, 3) + 10.* (x(:, 7) + x(:, 8));
    g = (g_temp+400.0551) +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1 - (y(:,1))./g);
  
   %constraints in objective space
    c(:,1) =max( -(y(:,1).^2 + y(:,2).^2-1), 0);
    c(:,2) = max(-( abs( (-y(:,1) + y(:,2) -0.5)/sqrt(2)) - 0.1/sqrt(2)), 0);
    c(:,3) = max(-( abs( (-y(:,1) + y(:,2) -0)/sqrt(2)) - 0.1/sqrt(2)), 0);
    c(:,4) = max(-( abs( (-y(:,1) + y(:,2) +0.5)/sqrt(2)) - 0.1/sqrt(2)), 0);
    
   %constraints in decision space
    c(:, 5) = x(:, 10).* x(:, 4) + 0.02.* x(:, 7) - 0.025.* x(:, 6);
    c(:, 6) = x(:, 10).* x(:, 5) + 0.02.* x(:, 8) - 0.015.* x(:, 9);
    c(:, 7) = abs(x(:, 2) + x(:, 3) - x(:, 4) - x(:, 5)) - 0.0001;
    c(:, 8) = abs(0.03.* x(:, 2) + 0.01.* x(:, 3) - x(:, 10).* (x(:, 4) + x(:, 5))) - 0.0001;
    c(:, 9) = abs(x(:, 4) + x(:, 7) - x(:, 6)) - 0.0001;
    c(:, 10) = abs(x(:, 5) + x(:, 8) - x(:, 9)) - 0.0001;
    y = y';
    c= -c';
end


function [y,c] = L4(x)


     x = x';
    % basic multi-objective problem
    g_temp = (x(:, 2) - 10).^2 + 5 * (x(:, 3) - 12).^2 + x(:, 4).^4 + 3 * (x(:, 5) - 11).^2 + 10 * x(:, 6).^6 + ...
            7 * x(:, 7).^2 + x(:, 8).^4 - 4 * x(:, 7).* x(:, 8) - 10 * x(:, 7) - 8 * x(:, 8);

    g = g_temp-680.6300573745 +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1-sqrt(y(:,1))./g);
  
   %constraints in objective space
    c(:,1) = max( -(y(:,1) + y(:,2)-1), 0);
    c(:,2) = max(- ( y(:,1)+ y(:,2) - 1 - abs(sin(10*pi*(y(:,1) - y(:,2) + 1) ))), 0);
    
    
   %constraints in decision space
    c(:, 3) = -127 + 2 * x(:, 2).^2 + 3 * x(:, 3).^4 + x(:, 4) + 4 * x(:, 5).^2 + 5 * x(:, 6);
    c(:, 4) = -282 + 7 * x(:, 2) + 3 * x(:, 3) + 10 * x(:, 4).^2 + x(:, 5) - x(:, 6);
    c(:, 5) = -196 + 23 * x(:, 2) + x(:, 3).^2 + 6 * x(:, 7).^2 - 8 * x(:, 8);
    c(:, 6) = 4 * x(:, 2).^2 + x(:, 3).^2 - 3 * x(:, 2).* x(:, 3) + 2 * x(:, 4).^2 + 5 * x(:, 7) - 11 * x(:, 8);
     
    y = y';
    c= -c';
end


function [y,c] = L5(x)

    x = x';
   % basic multi-objective problem
    g_temp = x(:, 2);
    g = g_temp-193.724510070035 +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1-sqrt(y(:,1))./g);
  
   %constraints in objective space
   
    c(:,1) = max( -(y(:,1) + y(:,2)-1), 0);
    c(:,2) = max(- ( y(:,1)+ y(:,2) - 1 - abs(sin(10*pi*(y(:,1) - y(:,2) + 1) ))), 0);
    c(:,3) = max( (y(:,1) - 0.8).*(y(:,2) - 0.6), 0); 
 
   %constraints in decision space
    c(:, 4) = -x(:, 2) + 35 * x(:, 3).^0.6 + 35 * x(:, 4).^0.6;
    c(:, 5) = abs(-300 * x(:, 4) + 7500 * x(:, 6) - 7500 * x(:, 7) - 25 * x(:, 5).* x(:, 6) + 25 * x(:, 5).* x(:, 7) + x(:, 4).* x(:, 5)) - 0.0001;
    c(:, 6) = abs(100 * x(:, 3) + 155.365 * x(:, 5) + 2500 * x(:, 8) - x(:, 3).* x(:, 5) - 25 * x(:, 5).* x(:, 8) - 15536.5) - 0.0001;
    c(:, 7) = abs(-x(:, 6) + log( - x(:, 5) + 900)) - 0.0001;
    c(:, 8) = abs(-x(:, 7) + log(x(:, 5) + 300)) - 0.0001;
    c(:, 9) = abs(-x(:, 8) + log(-2 * x(:, 5) + 700)) - 0.0001;
    y = y';
    c= -c';
end


function [y,c] = L6(x)
   
     x = x';
    % basic multi-objective problem
    g_temp = x(:, 2).^2 + x(:, 3).^2 + x(:, 2).* x(:, 3) - 14 * x(:, 2) - 16 * x(:, 3) + (x(:, 4) - 10).^2 + 4 * (x(:, 5) - 5).^2 + ...
            (x(:, 6) - 3).^2 + 2 * (x(:, 7) - 1).^2 + 5 * x(:, 8).^2 + 7 * (x(:, 9) - 11).^2 + 2 * (x(:, 10) - 10).^2 + (x(:, 11) - 7).^2 + 45;
    g = g_temp - 24.3062090681 +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1-sqrt(y(:,1))./g);
  
   %constraints in objective space
    c(:,1) = max( -(y(:,1) + y(:,2)-1), 0);
    c(:,2) = max(-(y(:,1) - 0.5).*( y(:,1)+ y(:,2) - 1 - abs(sin(10*pi*(y(:,1) - y(:,2) + 1) ))  ), 0);
    
   %constraints in decision space
    c(:, 3) = -105 + 4 * x(:, 2) + 5 * x(:, 3) - 3 * x(:, 8) + 9 * x(:, 9);
    c(:, 4) = 10 * x(:, 2) - 8 * x(:, 3) - 17 * x(:, 8) + 2 * x(:, 9);
    c(:, 5) = -8 * x(:, 2) + 2 * x(:, 3) + 5 * x(:, 10) - 2 * x(:, 11) - 12;
    c(:, 6) = 3 * (x(:, 2) - 2).^2 + 4 * (x(:, 3) - 3).^2 + 2 * x(:, 4).^2 - 7 * x(:, 5) - 120;
    c(:, 7) = 5 * x(:, 2).^2 + 8 * x(:, 3) + (x(:, 4) - 6).^2 - 2 * x(:, 5) - 40;
    c(:, 8) = x(:, 2).^2 + 2 * (x(:, 3) - 2).^2 - 2 * x(:, 2).* x(:, 3) + 14 * x(:, 6) - 6 * x(:, 7);
    c(:, 9) = 0.5 * (x(:, 2) - 8).^2 + 2 * (x(:, 3) - 4).^2 + 3 * x(:, 6).^2 - x(:, 7) - 30;
    c(:, 10) = -3 * x(:, 2) + 6 * x(:, 3) + 12 * (x(:, 10) - 8).^2 - 7 * x(:, 11);
    y = y';
    c= -c';
end

function [y,c] = L7(x)
      popsize = size(x,2);
     x = x';
    c1 = [-6.089 -17.164 -34.054 -5.914 -24.721 -14.986 -24.1 -10.708 -26.662 -22.179];
    % basic multi-objective problem
    g_temp = sum(x(:,2:11).* (repmat(c1, popsize, 1) + log(1E-30 + x(:,2:11)./repmat(1E-30 + sum(x(:,2:11), 2), 1, 10))), 2);
    g = g_temp +47.7648884595 +1;
    y(:,1) = x(:,1);
    y(:,2) = g.*(1-sqrt(y(:,1))./g);
  
   %constraints in objective space
    c(:,1) = max( -(y(:,1) + y(:,2)-1), 0);
    c(:,2) = max(-(y(:,1) - 0.5).*( y(:,1)+ y(:,2) - 1 - abs(sin(10*pi*(y(:,1) - y(:,2) + 1) ))  ), 0);
    c(:,3) = max(- ( abs(- y(:,1) + y(:,2))./sqrt(2) - 0.1./sqrt(2)), 0);
    
   %constraints in decision space
    c(:,4) = abs(x(:, 2) + 2 * x(:, 3) + 2 * x(:, 4) + x(:, 7) + x(:, 11) - 2) - 0.0001;
    c(:,5) = abs(x(:, 5) + 2 * x(:, 6) + x(:, 7) + x(:, 8) - 1) - 0.0001;
    c(:,6) = abs(x(:, 4) + x(:, 8) + x(:, 9) + 2 * x(:, 10) + x(:, 11) - 1) - 0.0001;
    y = y';
    c= -c';
end


function [y,c] = L8(x)

     x = x';
    % basic multi-objective problem
    g_temp = x(:, 3) + x(:, 4) + x(:, 5);
    g = g_temp-7049.2480205286 +1;
    
    y(:,1) = (x(:,1).*x(:,2)).*g;
    y(:,2) = (x(:,1).*(1 - x(:,2))).*g;
    y(:,3) = (1-x(:,1)).*g;
    
   %constraints in objective space
    c(:,1) = max( - (y(:,3) - 0.4).*(y(:,3) - 0.6), 0);
    
   %constraints in decision space
    c(:, 2) = -1 + 0.0025 * (x(:, 6) + x(:, 8));
    c(:, 3) = -1 + 0.0025 * (x(:, 7) + x(:, 9) - x(:, 6));
    c(:, 4) = -1 + 0.01 * (x(:, 10) - x(:, 7));
    c(:, 5) = -x(:, 3).* x(:, 8) + 833.33252 * x(:, 6) + 100 * x(:, 3) - 83333.333;
    c(:, 6) = -x(:, 4).* x(:, 9) + 1250 * x(:, 7) + x(:, 4).* x(:, 6) - 1250 * x(:, 6);
    c(:, 7) = -x(:, 5).* x(:, 10) + 1250000 + x(:, 5).* x(:, 7) - 2500 * x(:, 7);
    y = y';
    c= -c';
end


function [y,c] = L9(x)
 
     x = x';
    % basic multi-objective problem
    g_temp =  -0.5 * (x(:, 3).* x(:, 6) - x(:, 4).* x(:, 5) + x(:, 5).* x(:, 11) - x(:, 7).* x(:, 11) + x(:, 7).* x(:, 10) - x(:, 8).* x(:, 9));
    g = g_temp +0.8660254038 +1;
    
    y(:,1) = cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2)).*g;  
    y(:,2) = cos(0.5*pi*x(:,1)).*sin(0.5*pi*x(:,2)).*g;
    y(:,3) = sin(0.5*pi*x(:,1)).*g;
   
   %constraints in objective space
     c(:, 1) = max( -( y(:,1).^2 + y(:,2).^2 - 1), 0);
    
   %constraints in decision space
     c(:, 2) = x(:, 5).^2 + x(:, 6).^2 - 1;
     c(:, 3) = x(:, 11).^2 - 1;
     c(:, 4) = x(:, 7).^2 + x(:, 8).^2 - 1;
     c(:, 5) = x(:, 3).^2 + (x(:, 4) - x(:, 11)).^2 - 1;
     c(:, 6) = (x(:, 3) - x(:, 7)).^2 + (x(:, 4) - x(:, 8)).^2 - 1;
     c(:, 7) = (x(:, 3) - x(:, 9)).^2 + (x(:, 4) - x(:, 10)).^2 - 1;
     c(:, 8) = (x(:, 5) - x(:, 7)).^2 + (x(:, 6) - x(:, 8)).^2 - 1;
     c(:, 9) = (x(:, 5) - x(:, 9)).^2 + (x(:, 6) - x(:, 10)).^2 - 1;
     c(:, 10) = x(:, 9).^2 + (x(:, 10) - x(:, 11)).^2 - 1;
     c(:, 11) = x(:, 4).* x(:, 5) - x(:, 3).* x(:, 6);
     c(:, 12) = -x(:, 5).* x(:, 11);
     c(:, 13) = x(:, 7).* x(:, 11);
     c(:, 14) = x(:, 8).* x(:, 9) - x(:, 7).* x(:, 10);
     y = y';
    c= -c';
end







