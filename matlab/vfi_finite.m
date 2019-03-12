function [V,Cstar,par] = vfi_finite(par)
    
    % 1. allocate memory
    V = cell(par.T,1);
    Cstar = cell(par.T,1);
    if par.rho < 1.0
        par.grid_M = nonlinspace(0,par.M_max,par.NM,1.2); 
    else
        par.grid_M = nonlinspace(1e-4,par.M_max,par.NM,1.2); 
    end
    
    % 2. last period
    [V{par.T},Cstar{par.T}] = find_V(par,1);
    
    % 3. backwards over time
    for t = par.T-1:-1:1
        
        % a. interpolant
        par.V_plus_interp = griddedInterpolant(par.grid_M,V{t+1},'linear');   
        
        % b. find V for all states
        [V{t},Cstar{t}] = find_V(par,0);             
        
    end
        
end