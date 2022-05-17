import neat
xor_inputs=[(0.0,0.0),(0.0,1.0),(1.0,0.0),(1.0,1.0)]
xor_outputs=[(0.0,),(1.0,),(1.0,),(0.0,)]
def eval_fitness(net):
    error_sum=0.0
    for xi,xo in zip(xor_inputs,xor_outputs):
        output=net.activate(xi)
        error_sum+=abs(output[0]-xo[0])
    fitness=(4-error_sum)**2
    return fitness





