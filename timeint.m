function [parm, flow] = timeint(parm, flow)

    % timintegration follows here
    
    % loop over timesteps 
    for itst = 1 : parm.ntst

        disp(['time step ', num2str(itst)]);

        % compute right hand side of Navier-Stokes
        % (note: pressure gradient term is practically neglected)
        [flow.rhsu, flow.rhsv] = rhs_ns(parm, flow);

        % compute provisional velocity field 
        % (note: flow.u and flow.v are updated)
        [flow] = runge_kutta_2d_vec(parm, flow);

        % compute pressure difference (dp) field based on the
        % provisional velocity
        [flow] = direct_press_corr(parm, flow);

        % finally correct the velocity field by projection
        % (note: pressure field is also updated)
        [flow] = project(parm, flow);

        % check if the updated velocity fields are indeed
        % divergence-free
        [flow] = div_cal(parm, flow);
        fprintf('DEBUG: max div of final vel. is %e\n', div_check(flow));

    end    
end

