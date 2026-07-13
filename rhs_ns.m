function [rhsu, rhsv] = rhs_ns(parm, flow)

    % computes rhs for cases CHANNEL and SHEAR
    
    % local initialization
    
    rhsu = zeros(parm.m, parm.n);
    rhsv = zeros(parm.m, parm.n);
    
    % computation of right hand side in a loop
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            % diffusion_u() & diffusion_v() 
            % can be found in a separate file
            
            ad = advection_u(parm, flow, i, j);
            df = diffusion_u(parm, flow, i, j);
            pg = pressgrad_x(parm, flow, i, j);
            gr = gravity_x(parm);
            
            rhsu(i,j) = -ad-(1/parm.rho)*pg+gr+df;
            
            ad = advection_v(parm, flow, i, j);
            df = diffusion_v(parm, flow, i, j);
            pg = pressgrad_y(parm, flow, i, j);
            gr = gravity_y(parm);
            
            rhsv(i,j) = -ad-(1/parm.rho)*pg+gr+df;
            
        end
    end
    
end

% used functions

function [adv_u] = advection_u(parm, flow, i, j)

    % returns the advective term for v at {x(i),y(j)}
    
    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    
    adv_u = flow.u(i,j)*(flow.u(ip,j)-flow.u(im,j))/(2*parm.dx) + ...
            flow.v(i,j)*(flow.u(i,jp)-flow.u(i,jm))/(2*parm.dy);

end

function [adv_v] = advection_v(parm, flow, i, j)

    % returns the advective term for v at {x(i),y(j)}
    
    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    
    adv_v = flow.u(i,j)*(flow.v(ip,j)-flow.v(im,j))/(2*parm.dx) + ...
            flow.v(i,j)*(flow.v(i,jp)-flow.v(i,jm))/(2*parm.dy);
        
end

function [dpdx] = pressgrad_x(parm, flow, i, j)

    % returns the pressure gradient dp/dx at {x(i),y(j)}

    ip = parm.ip(i);
    im = parm.im(i);
    
    dpdx = (flow.p(ip,j) - flow.p(im,j))/(2*parm.dx);
    

end

function [dpdy] = pressgrad_y(parm, flow, i, j)

    % returns the pressure gradient dp/dy at {x(i),y(j)}

    jp = parm.jp(j);
    jm = parm.jm(j);
    
    dpdy = (flow.p(i,jp) - flow.p(i,jm))/(2*parm.dy);

end

function [gx] = gravity_x(parm)

    gx = parm.grav * (1);

end

function [gy] = gravity_y(parm)

    gy = parm.grav * 0;

end









