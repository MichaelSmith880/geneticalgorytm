import matplotlib.pyplot as plt
import random
import math


# Calculated function
def f(args):
    return f2(args)


def f1(args):
    return (3 - (math.sin(2 * args[0])) ** 2 - (math.sin(2 * args[1])) ** 2)


def f2(args):
    x = 1
    for i in range(len(args)):
        z = 0
        for j in range(5):
            z += (j + 1) * math.cos(((j + 1) + 1) * args[i] + (j + 1))
        x *= z
    return x


# Adaptation function
def s(x):
    return s2(x)


def s1(x):
    return math.exp(-abs(x - 1))


def s2(x):
    return math.exp(-abs(x + 187))


# Calculate the value of the 2-turned sequence represented by
'''
 Decoding and costing
 Chromosome group
 Chromosome Chrom_legth Length
 Max_Value, Min_Value upper and lower limit
 Div borders
'''


def b2d(b, chrom_length, max_value, min_value, div):
    rwno = []
    # Because there are several variables in the chromosome, it is necessary to separate
    for i in range(len(div)):
        if i == 0:
            star = 0
            end = div[i]
        else:
            star = div[i - 1] + 1
            end = div[i]
        t = 0
        for j in range(star, end):  # [[1, 2, 3 || 4, 5, 6]
            t += b[j] * (math.pow(2, j - star))
        t = t * max_value / (math.pow(2, end - star + 1) - 1) - min_value
        rwno.append(t)

    return rwno  # This is the list


'''
 Calculate the current value of a function
 Chromosome group
 Chromosome Chrom_legth Length
 Max_value, min_value minimum minimum
 Various segmentation
'''


def calobjValue(group, chrom_length, max_value, min_value, divid):
    obj_value = []
    for i in range(len(group)):
        x = b2d(group[i], chrom_length, max_value, min_value, divid)  #It can be multiple variables
        obj_value.append(f(x))
    return obj_value


# Get adaptation value
def calfitValue(obj_value):
    fit_value = []
    for i in range(len(obj_value)):
        temp = s(obj_value[i])  # Call the adaptation function calculation
        fit_value.append(temp)
    return fit_value


#The cumulative adaptation value is useful for calculating the average
def sum_fit(fit_value):
    total = 0
    for i in range(len(fit_value)):
        total += fit_value[i]
    return total


#
def selection(group, fit_value):
    newfit_value = []  # [[[[Chromosome], [anchor]], ...]
    newgroup = []  # []], [mother], [parent], [mother], ....]
    # Getage sum
    total_fit = sum_fit(fit_value)
    # Set each anchor point
    t = 0
    for i in range(len(group)):
        t += fit_value[i] / total_fit
        newfit_value.append([group[i], t])
    #            
    for i in range(len(newfit_value)):
        parents = len(newfit_value)  # initialization pointer
        r = random.random()  # указатель
        for j in range(len(newfit_value)):  #
            if newfit_value[j][1] > r:
                parents = j
                break
        newgroup.append(newfit_value[parents][0])

    return newgroup


#    .
def crossover(group, fit_value, pc):
    parents_group = selection(group, fit_value)  # [ [[parents]],....]
    group_len = len(parents_group)
    for i in range(0, group_len, 2):
        if (random.random() < pc):  # See if you want to compose
            cpoint = random.randint(0, len(parents_group[0]))  # random intersection
            temp1 = []
            temp2 = []
            temp1.extend(parents_group[i][0:cpoint])
            temp1.extend(parents_group[i + 1][cpoint:len(parents_group[i])])
            temp2.extend(parents_group[i + 1][0:cpoint])
            temp2.extend(parents_group[i][cpoint:len(parents_group[i])])
            group[i] = temp1
            group[i + 1] = temp2


# mutation gene
def mutation(group, pm):
    px = len(group)
    py = len(group[0])

    for i in range(px):  # intersection
        if (random.random() < pm):
            mpoint = random.randint(0, py - 1)  #
            if (group[i][mpoint] == 1):
                group[i][mpoint] = 0
            else:
                group[i][mpoint] = 1


'''
 Find the Best Solution and Optimal Solution Gene Making
 Group population
 Fit_Value Population Adaptation
'''


def best(group, fit_value):
    px = len(group)
    best_in = group[0]
    best_fit = fit_value[0]
    for i in range(1, px):
        if (fit_value[i] > best_fit):
            best_fit = fit_value[i]
            best_in = group[i]
    # print(best_in)
    return [best_in, best_fit]


'''
 Create initial population
 Population size Group_size
 Chromosome Chrom_legth Length
'''


def getFisrtGroup(group_size, chrom_length):
    # print ("initial population:")
    group = []
    for i in range(group_size):
        temp = []
        for j in range(chrom_length):
            temp.append(random.randint(0, 1))
        group.append(temp)
    # print(group)

    return group


generation = 50  # (less number, from results, less iteration)
group_size = 400  # , even
max_value = 20  # Sphere
min_value = 10  # Offset Correction
chrom_length = 800  #
divid = [399, chrom_length - 1]  # Enter value split point, last bit must be chromosome length
pc = 0.7  #
pm = 0.1  # variational probability
results = []  # Keep the best solution for every generation
fit_value = []  # Individual fitness
points = []  # More optimal solutions
#      
group = getFisrtGroup(group_size, chrom_length)

for i in range(generation):
    if i > 100:
        pm = 0.01
    if i > 1000:
        pm = 0.001
    obj_value = calobjValue(group, chrom_length, max_value, min_value, divid)  # Individual assessment
    fit_value = calfitValue(obj_value)  # Get group adaptation value
    best_individual, best_fit = best(group, fit_value)  # Return to optimal gene, optimal adaptation value

    xx = b2d(best_individual, chrom_length, max_value, min_value, divid)
    if (abs(f(xx) + 186.730909) < 0.000001):  # Find the best solution
        flag = False
        for p in points:
            if ((abs(xx[0] - p[0]) < 0.1) and (abs(xx[1] - p[1]) < 0.1)):  #
                flag = True
                break
        if flag == False:
            print(xx)
            points.append(xx)

    results.append([i, best_fit, b2d(best_individual, chrom_length, max_value, min_value, divid), best_individual])  #
    crossover(group, fit_value, pc)  # .
    mutation(group, pm)  # Mutations

# results.sort(key=lambda x:x[1])

rank = sorted(results, key=lambda x: x[1])
# print('\n', rank[-1])

# print(results)
x = b2d(rank[-1][3], chrom_length, max_value, min_value, divid)
# Final result
print("f(x) = ", f(x), "x = ", x, "Chromosome =", rank[-1][3], "Importance of adaptation =", rank[-1][1], "Algebra: ", rank[-1][0])

#              
X = []
Y = []
for i in range(generation):
    X.append(i)
    Y.append(results[i][1])

plt.plot(X, Y)

plt.show()
