import matplotlib.pyplot as plt
import random
import math


# Рассчитанная функция
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


# Функция адаптации
def s(x):
    return s2(x)


def s1(x):
    return math.exp(-abs(x - 1))


def s2(x):
    return math.exp(-abs(x + 187))


# Рассчитайте значение 2-овернутой последовательности, представленной
'''
 Декодирование и расчет стоимости
 Группа хромосома
 Хромосома Chrom_legth Длина
 Max_Value, Min_Value верхний и нижний предел
 Границы Div
'''


def b2d(b, chrom_length, max_value, min_value, div):
    rwno = []
    # Потому что в хромосоме есть несколько переменных, необходимо разделить
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

    return rwno  # Это список


'''
 Рассчитайте текущее значение функции
 Группа хромосома
 Хромосома Chrom_legth Длина
 Max_value, min_value minimum minimum
 Различная сегментация
'''


def calobjValue(group, chrom_length, max_value, min_value, divid):
    obj_value = []
    for i in range(len(group)):
        x = b2d(group[i], chrom_length, max_value, min_value, divid)  # Это может быть несколько переменных

        obj_value.append(f(x))
    return obj_value


# Получить значение адаптации
def calfitValue(obj_value):
    fit_value = []
    for i in range(len(obj_value)):
        temp = s(obj_value[i])  # Позвоните в расчет функции адаптации
        fit_value.append(temp)
    return fit_value


# Совокупное значение адаптации удобно для расчета среднего
def sum_fit(fit_value):
    total = 0
    for i in range(len(fit_value)):
        total += fit_value[i]
    return total


#
def selection(group, fit_value):
    newfit_value = []  # [[[[Хромосома], [якорь]], ...]
    newgroup = []  # []], [Мать], [родитель], [мать], ....]
    # Getage сумма
    total_fit = sum_fit(fit_value)
    # Установите каждую точку привязки
    t = 0
    for i in range(len(group)):
        t += fit_value[i] / total_fit
        newfit_value.append([group[i], t])
    #            
    for i in range(len(newfit_value)):
        parents = len(newfit_value)  # Указатель инициализации
        r = random.random()  # указатель
        for j in range(len(newfit_value)):  #
            if newfit_value[j][1] > r:
                parents = j
                break
        newgroup.append(newfit_value[parents][0])

    return newgroup


#    .
def crossover(group, fit_value, pc):
    parents_group = selection(group, fit_value)  # [ [[родители]],....]
    group_len = len(parents_group)
    for i in range(0, group_len, 2):
        if (random.random() < pc):  # Посмотреть, если вы хотите составить
            cpoint = random.randint(0, len(parents_group[0]))  # Случайное пересечение
            temp1 = []
            temp2 = []
            temp1.extend(parents_group[i][0:cpoint])
            temp1.extend(parents_group[i + 1][cpoint:len(parents_group[i])])
            temp2.extend(parents_group[i + 1][0:cpoint])
            temp2.extend(parents_group[i][cpoint:len(parents_group[i])])
            group[i] = temp1
            group[i + 1] = temp2


# Ген мутации
def mutation(group, pm):
    px = len(group)
    py = len(group[0])

    for i in range(px):  # Пересекание
        if (random.random() < pm):
            mpoint = random.randint(0, py - 1)  #
            if (group[i][mpoint] == 1):
                group[i][mpoint] = 0
            else:
                group[i][mpoint] = 1


'''
 Найти лучшее решение и оптимальное решение Gene решений
 Групповое население
 Адаптация населения Fit_Value
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
 Создать начальную популяцию
 Размер популяции Group_size
 Хромосома Chrom_legth Длина
'''


def getFisrtGroup(group_size, chrom_length):
    # print («начальное население:»)
    group = []
    for i in range(group_size):
        temp = []
        for j in range(chrom_length):
            temp.append(random.randint(0, 1))
        group.append(temp)
    # print(group)

    return group


generation = 50  # (меньшее количество, из результатов, менее итерация)
group_size = 400  # , даже
max_value = 20  # Сфера
min_value = 10  # Коррекция смещения
chrom_length = 800  #
divid = [399, chrom_length - 1]  # Введите точку разделения значения, последний бит должен быть длиной хромосомой
pc = 0.7  #
pm = 0.1  # Вариационная вероятность
results = []  # Храните лучшее решение для каждого поколения
fit_value = []  # Индивидуальная фитнес
points = []  # Более оптимальные решения
#      
group = getFisrtGroup(group_size, chrom_length)

for i in range(generation):
    if i > 100:
        pm = 0.01
    if i > 1000:
        pm = 0.001
    obj_value = calobjValue(group, chrom_length, max_value, min_value, divid)  # Индивидуальная оценка
    fit_value = calfitValue(obj_value)  # Получить ценность адаптации группы
    best_individual, best_fit = best(group, fit_value)  # Вернуться к оптимальному гену, оптимальное значение адаптации

    xx = b2d(best_individual, chrom_length, max_value, min_value, divid)
    if (abs(f(xx) + 186.730909) < 0.000001):  # Найти лучшее решение
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
    mutation(group, pm)  # Мутации

# results.sort(key=lambda x:x[1])

rank = sorted(results, key=lambda x: x[1])
# print('\n', rank[-1])

# print(results)
x = b2d(rank[-1][3], chrom_length, max_value, min_value, divid)
# Конечный результат
print("f(x) = ", f(x), "x = ", x, "Хромосома =", rank[-1][3], "Значение адаптации =", rank[-1][1], "Алгебра: ", rank[-1][0])

#              
X = []
Y = []
for i in range(generation):
    X.append(i)
    Y.append(results[i][1])

plt.plot(X, Y)

plt.show()