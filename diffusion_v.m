function [diff_v] = diffusion_v(parm, flow, i, j)

    % returns the diffusive term for v at {x(i),y(j)}, i.e. the Laplacian
    % multiplied with the viscosity

    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    dx = parm.dx;
    dy = parm.dy;
    nu = parm.nu;
    
    diff_v = nu*((flow.v(im,j) - 2*flow.v(i,j) + flow.v(ip,j))/dx^2 + ...
                 (flow.v(i,jm) - 2*flow.v(i,j) + flow.v(i,jp))/dy^2);
             
end