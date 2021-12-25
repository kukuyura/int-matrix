# Функция, которая ищет знак перестановки
def sign(perm):
    inv, n = 0, len(perm)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inv += 1
    return (-1) ** inv


# Функция, генерирующая перестановки списка
def perms(a):
    if len(a) == 1:
        return [a]
    ret = []
    for i in range(len(a)):
        for x in perms(a[:i] + a[i + 1:]):
            ret.append([a[i]] + x)
    return ret


# Функция, генерирующая сочетания размера k из списка (k <= len(a))
def combs(a, k):
    if k == 0:
        return []
    if k == 1:
        return [[x] for x in a]
    ret = []
    for i in range(len(a) - k + 1):
        for x in combs(a[i + 1:], k - 1):
            ret.append([a[i]] + x)
    return ret


# Класс матрицы, умющей работать только с целыми числами
class IntMatrix:
    # Матрица задается непустым двумерным массивом целых чисел
    def __init__(self, arr):
        self.a = arr
        self.n, self.m = len(arr), len(arr[0])

    # Складывать можно только матрицы одного размера
    def __add__(self, other):
        ret = [[0 for j in range(self.m)] for i in range(self.n)]
        for i in range(self.n):
            for j in range(self.m):
                ret[i][j] = self.a[i][j] + other.a[i][j]
        return IntMatrix(ret)

    # Вычитать можно только матрицы одного размера
    def __sub__(self, other):
        ret = [[0 for j in range(self.m)] for i in range(self.n)]
        for i in range(self.n):
            for j in range(self.m):
                ret[i][j] = self.a[i][j] - other.a[i][j]
        return IntMatrix(ret)

    # Матрицу n*m можно умножать только на матрицу m*k, либо на целое число
    def __mul__(self, other):
        if type(other) == int:
            ret = [[0 for j in range(self.m)] for i in range(self.n)]
            for i in range(self.n):
                for j in range(self.m):
                    ret[i][j] = self.a[i][j] * other
            return IntMatrix(ret)
        elif type(other) == IntMatrix:
            ret = [[0 for j in range(other.m)] for i in range(self.n)]
            for i in range(self.n):
                for j in range(other.m):
                    for k in range(self.m):  # self.m == other.n
                        ret[i][j] += self.a[i][k] * other.a[k][j]
            return IntMatrix(ret)

    def __rmul__(self, other):
        if type(other) == int:
            ret = [[0 for j in range(self.m)] for i in range(self.n)]
            for i in range(self.n):
                for j in range(self.m):
                    ret[i][j] = self.a[i][j] * other
            return IntMatrix(ret)
        elif type(other) == IntMatrix:
            ret = [[0 for j in range(other.n)] for i in range(self.m)]
            for i in range(other.n):
                for j in range(self.m):
                    for k in range(other.m):  # self.m == other.n
                        ret[i][j] += other.a[i][k] * self.a[k][j]
            return IntMatrix(ret)

    # Транспонирование
    def t(self):
        return IntMatrix([[self.a[i][j]
                           for i in range(self.n)] for j in range(self.m)])

    # Определитель можно узнать только у квадратной матрицы
    def det(self):
        perm_list = list(perms(list(range(self.n))))
        ret = 0
        for curr_perm in perm_list:
            local_product = sign(curr_perm)
            for i in range(self.n):
                local_product *= self.a[i][curr_perm[i]]
            ret += local_product
        return ret

    # Минор, в котором удалены строки c индексами exclude_rows
    # и столбцы с индексами exclude_cols,
    # можно узнать только для квадратной матрицы
    def minor(self, exclude_rows, exclude_cols):
        k = len(exclude_rows)
        ret = []
        for i in range(self.n):
            if i in exclude_rows:
                continue
            ret.append([])
            for j in range(self.n):
                if j in exclude_cols:
                    continue
                ret[-1].append(self.a[i][j])
        return IntMatrix(ret)

    # Характеристический многочлен можно узнать только у квадратной матрицы
    def characteristic_poly_coefs(self):
        coefs = [0 for i in range(self.n + 1)]
        for k in range(self.n):
            comb_list = list(combs(list(range(self.n)), k))
            for curr_comb in comb_list:
                coefs[k] += self.minor(curr_comb, curr_comb).det()
            coefs[k] *= (-1) ** (self.n - k)
        coefs[self.n] = 1
        return coefs

    # Присоединенную матрицу можно узнать только у квадратной матрицы
    def adj(self):
        ret = [[0 for j in range(self.n)] for i in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                ret[j][i] = self.minor([i], [j]).det()
        return IntMatrix(ret)

    # Строковое представление матрицы
    def __str__(self):
        return '\n'.join([' '. join([str(y) for y in x]) for x in self.a])
