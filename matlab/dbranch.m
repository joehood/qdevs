
function [ d ] = dbranch(R, L, E, iij, vi, vj)
    d = 1/L * (E - iij * R - vi - vj);
end