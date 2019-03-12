import numpy as np
import matplotlib.pyplot as plt

# %%

def find_V(par, last):
    """
    par: (object)
        parameter struct
    last : (boolean)
        if it is last period (True), else (False)
    """
    V, Cstar = np.empty(par.NM), np.empty(par.NM)

    # loop over states:
    for i_M in range(par.grid_M):

        Mt = par.grid_M[i_M]

        # value-of-choice function
        if last is True: #last period
            Vfunc = lambda C: par.u(C, par)
        else:
            Vfunc = lambda C: par.u(C,par) + par.beta*sum(par.w*par.V_plus_interp(par.R*(Mt-C)+par.Y))

        Vfunc_neq = lambda C: -Vfunc(C)

        if i_M == 1:
            initial_guess = Mt / 2
        else:
            initial_guess = Cstar(i_M - 1)

        # Optimizer


        # Save optimum
        V(i_M) = #VALUE OF valuefunciton



# %%

def vfi_finite(par):

    V = np.empty(par.T)
    Cstar = np.empty(par.T)

    if par.rho < 1.0:
        par.grid_M = np.linspace(0, par.M_max, par.NM)
    else:
        par.grid_M = np.linspace(1e-4, par.M_max, par.NM)

    # last period
    V[par.T], Cstar[par.T] = find_V(par, 1)

    return V, Cstar, par
# %%
