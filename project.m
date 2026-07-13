function [flow] = project(parm, flow)
% Projection step for obtaining a divergence-free velocity field, 
% i.e. updates pressure and velocity by
%                                                                                                                                                                                                                                                                                         
%  p = p + dp
%
%                      d (dp)
%  u = u - dt / rho * --------
%                        dx
%
%                      d (dp)
%  v = v - dt / rho * --------
%                        dy
%
%
% Again, watch out for the boundaries and do not update u or v if 
% par->atbound[j] == 1.

    for i = 1 : parm.m
        
        ip = parm.ip(i);
        im = parm.im(i);
    
        for j = 1 : parm.n
            
            jp = parm.jp(j);
            jm = parm.jm(j);
            
            % projection step follows here
            
            % flow.p(i,j) = ???
            
            % flow.u(i,j) = ???
                        
            % flow.v(i,j) = ???   
            flow.p(i,j) = flow.p(i,j) + flow.dp(i,j);

            if parm.atbounds(j) ~= 1
                flow.u(i,j) = flow.u(i,j) - parm.dt / parm.rho * (flow.dp(ip,j) - flow.dp(im,j)) / (2*parm.dx);
                flow.v(i,j) = flow.v(i,j) - parm.dt / parm.rho * (flow.dp(i,jp) - flow.dp(i,jm)) / (2*parm.dy);
            end 
        end
    end
    
end

