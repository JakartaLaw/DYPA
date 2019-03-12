function [V,Cstar] = find_V(par,last)
    
    % loop over states    
    V     = nan(par.NM,1);
    Cstar = nan(par.NM,1);    
    for i_M = 1:numel(par.grid_M)
    
        Mt = par.grid_M(i_M);   
        
        % a. value-of-choice function
        if last == 1 % last period
            Vfunc = @(C) par.u(C,par);
        else
            Vfunc = @(C) par.u(C,par) + par.beta*sum(par.w.*par.V_plus_interp(par.R*(Mt-C)+par.Y));
        end
        Vfunc_neg = @(C) -Vfunc(C); % negative due to minimizer
        
        % b. initial guess
        if i_M == 1
            initial_guess = Mt/2;
        else
            initial_guess = Cstar(i_M-1);
        end
        
        % c. find optimum
        [x,fval] = fmincon(Vfunc_neg,initial_guess,[],[],[],[],...
                           0,Mt,[],par.options);                     
        
        % d. save optimum       
        V(i_M)     = -fval;
        Cstar(i_M) = x;
    
    end  
        
end