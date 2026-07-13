function [divmax] = div_check(flow)

    % computes maximum divergence in field
    
    divmax = max(abs(flow.div(:)));

end
