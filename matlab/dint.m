
function [ x ] = dint(t, tlast, x0, d)

    % update elapsed time for each component:
    dt = t - tlast;

    % update internal states:
    x = x0 + d * dt;
    
    
end

