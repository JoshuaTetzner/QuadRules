using Plots
pyplot()


order = 3
x, w = asymmetriccubature(order)
plot(x[:,1], x[:,2], w[:], st:=surface)

