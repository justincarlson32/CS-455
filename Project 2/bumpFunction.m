function bump = bumpFunction(z)
    if (z >= 0) && (z < .2) % 0->.2 results in bad smoothing
        bump = 1;
    elseif (z >= .2) && (z <= 1) % normal smoothing
        bump = .5*(1 + cos(pi*((z - .2) / .8)));
    else
        bump = 0;
    end
end