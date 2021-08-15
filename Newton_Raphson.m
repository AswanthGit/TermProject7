function sol = Newton_Raphson(q,n,tol,iter_max,fun,t)
% q = position coordinates at time t
% n = number of variables
% tol_max = maximum tolerance that is to be specified
% fun = function to be called for getting constraints and Jacobian Matrix
% t= current time step 
coords = zeros(iter_max,n+1);
flag = 0;sol = [];
for i = 1:iter_max
    [Phi,D] = fun(t,q);
    D = D(:,1:20);
    err = sqrt(Phi'*Phi);
    coords(i,:) = [err q'];
    if err<tol
       flag = 1;
       sol = coords(i,2:n+1)';
       break;
    end
    delta_q = -D\Phi;
    delta_q = [delta_q(1:20,1)' 0]';
    q = q + delta_q;
end
if flag == 0
    'Convergence failed in Newton-Raphson'
    return;
end
end