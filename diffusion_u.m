function [diff_u] = diffusion_u(parm, flow, i, j)

    % returns the diffusive term for u at {x(i),y(j)}, i.e. the Laplacian
    % multiplied with the viscosity
    
    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    dx = parm.dx;
    dy = parm.dy;
    nu = parm.nu;
    
    diff_u = nu*((flow.u(im,j) - 2*flow.u(i,j) + flow.u(ip,j))/(dx^2) + ...
                 (flow.u(i,jm) - 2*flow.u(i,j) + flow.u(i,jp))/(dy^2));      
             
end