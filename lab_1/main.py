from intervals import *

def det(delta, coef):
    A = [[Interval(coef[0] - delta, coef[0] + delta), Interval(coef[1] - delta, coef[1] + delta)],
      [Interval(coef[2] - delta, coef[2] + delta), Interval(coef[3] - delta, coef[3] + delta)]]
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def minDelta(a, b, eps, coef):
    while b - a > eps:
        x = a + (b - a) / 2
        detA = det(x, coef)
        if detA.a < 0 and detA.b > 0:
             b = x
        else:
            a = x
    detA = det(b, coef)
    return b, detA
    

if __name__ == "__main__":
    eps = 0.0000001

    coef = [float(n) for n in input().split()]

    delta, detA = minDelta(0, min(coef), eps, coef)

    print(f"delta = {delta}, ({detA.a}, {detA.b})")
    


        
