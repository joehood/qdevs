
function [ d ] = dnode(R, C, H, ii, ij, vi)
    d = 1/C * (H - vi / R + ii - ij);
end

