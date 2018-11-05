

function [ tnext ] = ta(d, x, q, dQ, t)
    if (d > 0)
        tnext = t + (q + dQ - x) / d;
                 
    elseif (d < 0)
        dt = (q - 0.5 * dQ - x) / d;
        tnext = t + abs(dt);
    else
        tnext = inf;
    end
    
     A=[d x q dQ t];
     
       if(tnext<0)
         disp(A);
%          display(d,x,q,dQ,t);
         error('tnext negavtive');
         
       end   
end

