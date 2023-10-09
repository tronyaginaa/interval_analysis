class Interval:

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __add__(self, other):
        return Interval(self.a + other.a, self.b + other.b)
    
    def __sub__(self, other):
        return Interval(self.a - other.b, self.b - other.a)
    
    def __mul__(self, other):
        return Interval(min(self.a * other.a, self.a * other.b, self.b * other.a, self.b * other.b), 
                        max(self.a * other.a, self.a * other.b, self.b * other.a, self.b * other.b))
    
    def __div__(self, other):
        return Interval(min(self.a / other.a, self.a / other.b, self.b / other.a, self.b / other.b), 
                        max(self.a / other.a, self.a / other.b, self.b / other.a, self.b / other.b))
