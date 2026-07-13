function [parm, flow] = initialize(parm, flow)

    % general initialisation: ip, im, jp, jm, atbounds 

    for i = 1 : parm.m

        parm.ip(i) = i+1;
        parm.im(i) = i-1;

    end
    
    for j = 1 : parm.n

        parm.jp(j) = j+1;
        parm.jm(j) = j-1;
        parm.atbounds(j) = 0;

    end
    
    % boundaries of domain (in the genral case)
    % specific definitions follow in the switch section
    
    parm.ip(parm.m) = 1;
    parm.jp(parm.n) = 1;
    parm.im(1) = parm.m;
    parm.jm(1) = parm.n;
    
    % dx, dy
    
    parm.dx = parm.xl/(parm.m-1);
    parm.dy = parm.yl/(parm.n-1);
    
    % set initial fields for different flow cases
    
    switch parm.flowtype
        
        case {'CHANNEL'}
            
            % gravity
            parm.grav = 9.81;
            
            % initial velocity field
            [parm, flow] = initcond_channel(parm, flow);
            
            % additional initialisations jm, jp, atbounds
            % for that specific case
            parm.jm(1) = 2;
            parm.jp(parm.n) = parm.n-1;
            parm.atbounds(1) = 1;
            parm.atbounds(parm.n) = 1;
            
            % Special initialisations of jpp, jmm for the Rhie-Chow interpolation
            % Only relevant for part 4b.
            % Ignore for part 4a!
            parm.jmm(1) = 3;
            parm.jpp(parm.n) = parm.n-2;
            parm.jmm(2) = 2;
            parm.jpp(parm.n-1) = parm.n-1;
            
            % initialize for Gauss Jordan
            [parm, flow] = init_gauss_jordan(parm, flow);
            
        case {'SHEAR'}
            
            % gravity
            parm.grav = 0;
            
            % initial velocity field
            [parm, flow] = initcond_shear(parm, flow);
            
            % initialize for Gauss Jordan
            [parm, flow] = init_gauss_jordan(parm, flow);
         
        otherwise 
            
            error('Flowtype not specified.')
    end

end

function [parm, flow] = initcond_channel(parm, flow)

% define initial velocity field u,v for the channel-case
% in the field and at boundaries

    for i = 1 : parm.m
        for j = 1 : parm.n
            
            % initialisation at field
            flow.u(i,j) = 0;
            % If you want to intensify the initial pertubation increase the
            % pertubation intensity prtb_int
            prtb_int = 0.01;
            flow.v(i,j) = prtb_int*sin(2.*pi*j/parm.n)*cos(pi*i/parm.m);
            flow.p(i,j) = 0;
            
        end
        
        % initialisation at boundary
        flow.u(i,1) = 0;
        flow.u(i,parm.n) = 0;
        flow.v(i,1) = 0;
        flow.v(i,parm.n) = 0;
        
    end
    
end

function [parm, flow] = initcond_shear(parm, flow)

% define initial velocity field u,v for the shear-case
% in the field (and at boundaries)

    for i = 1 : parm.m
        for j = 1 : parm.n
            
            % velocity component v, whole field
            flow.v(i,j) = 0.05*sin(getx(parm,i)*2*pi/parm.xl);
            flow.p(i,j) = 0;
            
        end
        
        for j = 1 : floor(parm.n/2)
                        
            % velocity component u, lower half
            flow.u(i,j) = tanh((4*gety(parm, j)-parm.yl)/parm.yl*7.5);
            
        end
        
        for j = ceil(parm.n/2) : parm.n
            
            % velocity component u, upper half
            flow.u(i,j) = tanh((3*parm.yl-4*gety(parm, j))/parm.yl*7.5);
        end
        
    end
    
end

function [parm, flow] = init_gauss_jordan(parm, flow)

    % initialize matrizes for Gauss-Jordan
    
    % if debugging is needed debugp = 1, else 0
    parm.debugp = 0;
    
    % copy some data locally for simplification   
    dt = parm.dt;
    rho = parm.rho;
    dx = parm.dx;
    dy = parm.dy;
    
    m = parm.m;
    n = parm.n;
    
    % Attention: BC's are important
    
    npoints = m*n;
    parm.npoints = npoints;

    % Matrix coefficients
    
    kx = (dt/rho) * (1/dx^2);
    kc = (dt/rho) * (2/dx^2 + 2/dy^2) * (-1);
    ky = (dt/rho) * (1/dy^2);
    
    % initialize coefficient matrix with zero
    
    flow.matrix_T  = zeros(parm.npoints, parm.npoints);

    % define matrix after initialization
    
    for i = 1 : m
        for j= 1 : n
            
            % row is the same for all elements
            ii = m*(j-1) + i;
            jj = ii;
            
            % center stencil element
            flow.matrix_T(ii,jj) = kc;
            
            % stencil elements in x-direction
            % im 
            jj = m*(j-1) + parm.im(i);
            flow.matrix_T(ii,jj) = kx;
            
            % ip
            jj = m*(j-1) + parm.ip(i);
            flow.matrix_T(ii,jj) = kx;
            
            % stencil elements in y-direction
            % jm
            jj = (m * (parm.jm(j)-1)) + i;
            flow.matrix_T(ii,jj) = flow.matrix_T(ii,jj) + ky;
            
            % jp
            jj = (m * (parm.jp(j)-1)) + i;
            flow.matrix_T(ii,jj) = flow.matrix_T(ii,jj) + ky;
               
        end
    end
            
    % de-singularize the system, by setting dp[1] = 0
    
    % when there is no explicit specification for a pressur level (e.g.
    % periodic BC), the pressure value can shift arbitrarily. Therefore we
    % set the value of the pressure in the first line so that the
    % system is determined! Best is zero ... then no values have to be
    % brought to the other side
    
    flow.matrix_T(1, 2:end) = 0;          % leads to dp[1] = 0
    
    % invert matrix for solving the linear system of equations 
    % in the function 'direct_press_corr.m'
    flow.matrix_invT = inv(flow.matrix_T);
    fprintf('matrix inversion is done\n')
    
end


%**************************************************************************

% helpful functions

function [x] = getx(parm, i)

% get x-value for a corresponding index i

    x = parm.xl/(parm.m-1)*(i-1);

end

function [y] = gety(parm, j)

% get y-value for a corresponding index j

    y = parm.yl/(parm.n-1)*(j-1);

end

