function [parm, flow] = set_params(parm, flow, infilename)

    load(infilename, 'flowtype', 'dt', 'ntst', 'm', 'n', ...
        'nu', 'rho', 'xl', 'yl');
    
    parm.flowtype = flowtype;
    parm.dt = dt;       
    parm.ntst = ntst;
    parm.m = m;         
    parm.n = n;
    parm.nu = nu;       
    parm.rho = rho;
    parm.xl = xl;       
    parm.yl = yl;
    
    %**********************************************************************
    % to manually override automatic input, uncomment below
    % and reset params:
    %
    % parm.dt = ...
    % parm.ntst = ...
    %
    %**********************************************************************
    
end




