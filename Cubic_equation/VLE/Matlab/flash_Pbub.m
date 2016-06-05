function flash_Pbub(T,init_guess)
%{
    Calculates bubble pressure(bar) at given temperauter(celcius)
    using init_guess to initiate iteration
    Assuming Ideal solution for binary component
    Add one Psat equation for any new component in the mixture
%} 

cmp1 = input('Enter component 1: ','s');
cmp2 = input('Enter Component 2: ','s');

% Get antoine coefficients
[cmp1,cmp2] = getAntione(cmp1,cmp2);

while(1)
    z = input('Enter feed composition in the order the name was given in a vector: ');
    if(length(z) > 1)
        break;
    end
end

% Component Psat equation
Psat1 = @(x) 10^(cmp1(1) - cmp1(2)/(x + cmp1(3)))/760*1.01325;
Psat2 = @(x) 10^(cmp2(1) - cmp2(2)/(x + cmp2(3)))/760*1.01325;

% iteration
opt = optimoptions('fsolve','display','iter','MaxIter',1e3);
%{
    Objective function:
    sum(yi) - 1 == 0
%}
func = @(x) z(1)*Psat1(T)/x + z(2)*Psat2(T)/x - 1; 
[P] = fzero(func,init_guess,opt); 
y1 = (z(1)*Psat1(T)/P);
y2 = (z(2)*Psat2(T)/P);

fprintf('Bubble Pressure = %.2f \n',P)
fprintf('y1 = %.3f | y2 = %.3f \n',y1,y2)


end

