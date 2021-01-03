from numpy import ones

# TODO. I can change these to zeros too to make plates set crust function easier


class HeightMap(ones):
    '''
        A class for constructing a matrix for height maps.
    '''

    def __init__(self, x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__name__)
        self.width = x
        self.height = y
        super().__init__((x, y), dtype=int)

    def set_all(self, v: int):
        for x in range(self.length):
            for y in range(x.length):
                self[x, y] = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height


class MassMap(ones):
    '''
        A class for constructing a matrix for mass maps.
    '''

    def __init__(self, x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__name__)
        self.width = x
        self.height = y
        super().__init__((x, y), dtype=float)

    def set_all(self, v: float):
        for x in range(self.length):
            for y in range(x.length):
                self[x, y] = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height


class AgeMap(ones):
    '''
        A class for constructing a matrix for age maps.
    '''

    def __init__(self,  x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__name__)
        self.width = x
        self.height = y
        super().__init__((x, y), dtype=int)

    def set_all(self, v: int):
        for x in range(self.length):
            for y in range(x.length):
                self[x, y] = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height
