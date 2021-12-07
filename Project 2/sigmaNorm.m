function sigmaNormFun = sigmaNorm(arg1)
    epsilon = 0.1; %(0,1)
    sigmaNormFun = (1/epsilon)*(sqrt(1 + epsilon*(norm(arg1,2)^2)) - 1);
end
