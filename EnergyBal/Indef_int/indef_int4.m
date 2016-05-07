    function [resp2, roots] = indef_int4(A ,string,lowerBound,HQ) 
    % Finds indefinite integral and calculates roots
    % returns symbols
    
       R = 8.314; %[=]J/mol/K
       LB = lowerBound; %
       syms T ;
       sym(HQ);
       
       % This Block Determines the right equation
        if(strcmpi(string,'CP'))
            
           eqtn = A;
           resp = int(eqtn); % Integral
           T = LB; % inittial temperature
           resp2 = resp - subs(resp) + HQ;
           roots = vpa(solve(resp2)); % gets roots
          
           
        elseif(strcmpi(string,'CV'))
                
            eqtn = A-R;        
            resp = int(eqtn);
            T = LB; % inittial temperature
            resp2 = resp - subs(resp) + HQ;
            roots = vpa(solve(resp2)); % gets roots