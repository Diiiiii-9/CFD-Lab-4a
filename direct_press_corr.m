function [flow] = direct_press_corr(parm, flow)

    % compute divergence of intermediate velocity field 
    % (note: flow.div is an [m,n]-matrix)
    [flow] = div_cal(parm, flow);
    fprintf('DEBUG: max div of int vel. is %e\n', div_check(flow));
    
    % reshaping divergence to [m*n,1]-1D array
    % this gives the rhs for the linear system of equations
    Mdiv = reshape(flow.div, [parm.npoints,1]);
    
    % pinning the pressure, because no pressure level has ben specified
    % compare comment in function initialize.m
    
    Mdiv(1) = 0;
    
    % solving Poisson-equation directly (linear system of equation, by matrix inversion)
    % flow.matrix_T * Mdp = Mdiv
    
    Mdp = flow.matrix_invT * Mdiv;
    
    % reshape Mdp 
    
    flow.dp = reshape(Mdp, parm.m, parm.n);
    
end

