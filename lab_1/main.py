from intervals import *

def det(delta, coef, rad):
    A = [[Interval(coef[0] - delta * rad[0], coef[0] + delta * rad[0]), Interval(coef[1] - delta* rad[1], coef[1] + delta* rad[1])],
      [Interval(coef[2] - delta* rad[2], coef[2] + delta* rad[2]), Interval(coef[3] - delta* rad[3], coef[3] + delta* rad[3])]]
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def minDelta(a, b, eps, coef, rad):
    while b - a > eps:
        x = a + (b - a) / 2
        detA = det(x, coef, rad)
        if detA.a < 0 and detA.b > 0:
             b = x
        else:
            a = x
    detA = det(b, coef, rad)
    return b, detA
    

if __name__ == "__main__":
    eps = 0.0000001

    print("Enter MAtrix A coefficients")
    coef = [float(n) for n in input().split()]
    print("Enter MAtrix M coefficients")
    rad = [float(n) for n in input().split()]

    delta, detA = minDelta(0, min(coef), eps, coef, rad)

    print(f"delta = {delta}, ({detA.a}, {detA.b})")
    


        
