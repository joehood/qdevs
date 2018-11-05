

function [ q ] = quantize(x, q, dQ, eps)

    if x >= q + dQ - eps
        q = q + dQ;
    elseif x <= q - 0.5 * dQ + eps
        q = q - dQ;
    end

end