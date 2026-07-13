function [flow] = runge_kutta_2d_vec(parm, flow)

    % 1st. RK step
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            flow.u(i,j) = flow.u(i,j) + (1/3)*parm.dt*flow.rhsu(i,j)*...
                (1-parm.atbounds(j));
            
            flow.v(i,j) = flow.v(i,j) + (1/3)*parm.dt*flow.rhsv(i,j)*...
                (1-parm.atbounds(j));
            
        end
    end
    
    % 2nd. RK step
    
    [flow.secrhsu, flow.secrhsv] = rhs_ns(parm, flow);
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            flow.secrhsu(i,j) = flow.secrhsu(i,j) -(5/9)*flow.rhsu(i,j);
            flow.secrhsv(i,j) = flow.secrhsv(i,j) -(5/9)*flow.rhsv(i,j);
            
            flow.u(i,j) = flow.u(i,j)+(15/16)*parm.dt*flow.secrhsu(i,j)*...
                (1-parm.atbounds(j));
            flow.v(i,j) = flow.v(i,j)+(15/16)*parm.dt*flow.secrhsv(i,j)*...
                (1-parm.atbounds(j));
        
        end
    end
    
    % 3rd. RK step
    
    [flow.rhsu, flow.rhsv] = rhs_ns(parm, flow);
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            flow.rhsu(i,j) = flow.rhsu(i,j)-(153/128)*flow.secrhsu(i,j);
            flow.rhsv(i,j) = flow.rhsv(i,j)-(153/128)*flow.secrhsv(i,j);
            
            flow.u(i,j) = flow.u(i,j)+(8/15)*parm.dt*flow.rhsu(i,j)*...
                (1-parm.atbounds(j));
            flow.v(i,j) = flow.v(i,j)+(8/15)*parm.dt*flow.rhsv(i,j)*...
                (1-parm.atbounds(j));
            
        end
    end
    
end