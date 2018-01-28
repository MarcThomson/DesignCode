%DiffEqSolver

initialConditions;
parameterDefinitions;

tspan = [1 24*7];

[t,y] = ode23s(@dydt,tspan,y0);

function derivative_system = dydt(t,y)
    V = 12000*1000; %mL
    extraCellularComponents;
    load('stoichMatrix.mat')
    %y = (y>=0).*y;
    
    inputs = mat2cell(y,ones(1,length(y)),1);
    rates = rateLaws(inputs {:});
    derivative_system = stoichMatrix * rates;
    
    for i = 1:length(extraCellular)
        if extraCellular(i)==1
            derivative_system(i) = derivative_system(i) * y(end) * 1000;  %mM/hr
        end
    end
    
end