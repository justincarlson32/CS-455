function point = pointFromAction(a)
    point = [0.,0.];
    if (a == 1)
        point = [540,400];
    elseif (a == 2)
        point = [540,0];
    elseif (a == 3)
        point = [300,400];
    elseif (a == 4)
        point = [300,0];
    end
end