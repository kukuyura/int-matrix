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


# Алгоритм Евклида
def gcd(a, b):
    if a < b:
        a, b = b, a
    while b != 0:
        a, b = b, a % b
    return a


# Класс матрицы, умеющей работать только с целыми числами
class IntMatrix:
    # Матрица задается непустым двумерным массивом целых чисел
    def __init__(self, arr):
        self.a = arr
        self.n, self.m = len(arr), len(arr[0])

    # Складывать можно только матрицы одного размера
    def __add__(self, other):
        ret = [[0 for _2 in range(self.m)] for _ in range(self.n)]
        for i in range(self.n):
            for j in range(self.m):
                ret[i][j] = self.a[i][j] + other.a[i][j]
        return IntMatrix(ret)

    # Вычитать можно только матрицы одного размера
    def __sub__(self, other):
        ret = [[0 for _2 in range(self.m)] for _ in range(self.n)]
        for i in range(self.n):
            for j in range(self.m):
                ret[i][j] = self.a[i][j] - other.a[i][j]
        return IntMatrix(ret)

    # Матрицу n*m можно умножать только на матрицу m*k, либо на целое число
    def __mul__(self, other):
        if type(other) == int:
            ret = [[0 for _2 in range(self.m)] for _ in range(self.n)]
            for i in range(self.n):
                for j in range(self.m):
                    ret[i][j] = self.a[i][j] * other
            return IntMatrix(ret)
        elif type(other) == IntMatrix:
            ret = [[0 for _2 in range(other.m)] for _ in range(self.n)]
            for i in range(self.n):
                for j in range(other.m):
                    for k in range(self.m):  # self.m == other.n
                        ret[i][j] += self.a[i][k] * other.a[k][j]
            return IntMatrix(ret)

    def __rmul__(self, other):
        if type(other) == int:
            ret = [[0 for _2 in range(self.m)] for _ in range(self.n)]
            for i in range(self.n):
                for j in range(self.m):
                    ret[i][j] = self.a[i][j] * other
            return IntMatrix(ret)
        elif type(other) == IntMatrix:
            ret = [[0 for _2 in range(other.n)] for _ in range(self.m)]
            for i in range(other.n):
                for j in range(self.m):
                    for k in range(other.m):  # self.m == other.n
                        ret[i][j] += other.a[i][k] * self.a[k][j]
            return IntMatrix(ret)

    # Функция бинарного возведения матрицы в положительную целую степень
    def __pow__(self, p):
        if p == 1:
            return self
        if p % 2:
            return self * self ** (p - 1)
        get = self ** (p // 2)
        return get * get

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

    # Минор, в котором удалены строки с индексами exclude_rows
    # и столбцы с индексами exclude_cols (нумерация с 0),
    # можно узнать только для квадратной матрицы
    def minor(self, exclude_rows, exclude_cols):
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
        coefs = [0 for _ in range(self.n + 1)]
        coefs[0] = (-1) ** self.n * self.det()
        for k in range(1, self.n):
            comb_list = list(combs(list(range(self.n)), k))
            for curr_comb in comb_list:
                coefs[k] += self.minor(curr_comb, curr_comb).det()
            coefs[k] *= (-1) ** (self.n - k)
        coefs[self.n] = 1
        return coefs

    # Присоединенную матрицу можно узнать только у квадратной матрицы
    def adj(self):
        ret = [[0 for _2 in range(self.n)] for _ in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                ret[j][i] = (-1) ** (i + j) * self.minor([i], [j]).det()
        return IntMatrix(ret)

    # Приведение к ступенчатому виду, может работать медленно на больших числах
    # x, y означают, что нужно запустить алгоритм с соответствующей строки и столбца
    # (нумерация с 0)
    def row_echelon_form(self, x=0, y=0):
        if x == self.n or y == self.m:
            return
        for i in range(x, self.n):
            if self.a[i][y] != 0:
                self.a[x], self.a[i] = [_ for _ in self.a[i]], \
                                       [_ for _ in self.a[x]]
                break
        else:
            self.row_echelon_form(x, y + 1)
            return
        for i in range(x + 1, self.n):
            coef = self.a[i][y]
            for j in range(y, self.m):
                self.a[i][j] = self.a[x][y] * self.a[i][j] - \
                               coef * self.a[x][j]
        self.row_echelon_form(x + 1, y + 1)
        for i in range(x):
            coef = self.a[i][y]
            for j in range(self.m):
                self.a[i][j] = self.a[x][y] * self.a[i][j] - \
                               coef * self.a[x][j]

    # Ранг матрицы, может работать медленно на больших числах
    def rk(self):
        copy = IntMatrix([[y for y in x] for x in self.a])
        copy.row_echelon_form()
        rk = self.n
        for i in range(copy.n):
            if copy.a[i] == [0 for _ in range(copy.m)]:
                rk -= 1
        return rk

    # Фундаментальная система решений - fundamental system of solutions
    # Возвращает ФСР, уложенную в столбцы матрицы, может работать медленно на больших числах,
    def fss(self):
        copy = IntMatrix([[y for y in x] for x in self.a])
        copy.row_echelon_form()
        main_vars_rows = []
        main_vars_cols = []
        coef = 1
        for i in range(copy.n):
            for j in range(copy.m):
                if copy.a[i][j] != 0:
                    main_vars_rows.append(i)
                    main_vars_cols.append(j)
                    coef *= copy.a[i][j]
                    break
        rk = len(main_vars_cols)
        ans = []
        for pos in range(copy.m):
            if pos in main_vars_cols:
                continue
            x = [0 for _ in range(copy.m)]
            x[pos] = coef
            for i in range(rk):
                x[main_vars_cols[i]] = -(copy.a[main_vars_rows[i]][pos] *
                                         coef // copy.a[main_vars_rows[i]][main_vars_cols[i]])
            ans.append(x)
        return IntMatrix(ans).t()

    # Делит каждую строку матрицы на ее НОД
    def divide_every_row_by_gcd(self):
        for i in range(self.n):
            d = 0
            for j in range(self.m):
                d = gcd(d, self.a[i][j])
            if d == 0:
                continue
            for j in range(self.m):
                self.a[i][j] //= d

    # x-ая строка матрицы (нумерация с 0)
    def get_row(self, x):
        return IntMatrix([[self.a[x][j] for j in range(self.m)]])

    # y-ый столбец матрицы (нумерация с 0)
    def get_col(self, y):
        return IntMatrix([[self.a[i][y]] for i in range(self.n)])

    # Строковое представление матрицы
    def __str__(self):
        return '[[' + '],\n ['.join([', '. join([str(y) for y in x]) for x in self.a]) + ']]'


# Единичная матрица размера n
def E(n):
    return IntMatrix([[1 if i == j else 0 for j in range(n)] for i in range(n)])
