function [flow] = div_cal(parm, flow)

    % approximates the divergence of the velocity field, using a 
    % central difference scheme
    % needed to solve the poisson equation
    % div = du/dx + dv/dy
    
    dx = parm.dx;
    dy = parm.dy;
    
    for i = 1 : parm.m
        
        ip = parm.ip(i);
        im = parm.im(i);
        
        for j = 1 : parm.n
            
            jp = parm.jp(j);
            jm = parm.jm(j);
            
            flow.div(i,j) = ((flow.u(ip,j)-flow.u(im,j))/(2*dx) + ...
                             (flow.v(i,jp)-flow.v(i,jm))/(2*dy)) * ...
                             (1 - parm.atbounds(j)); 
        
        end
    end
    
end