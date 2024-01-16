import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt

def create_L(X):
    midL = np.zeros((2, 2))
    radL = np.zeros((2, 2))
    radL[0][0] = 0
    midL[0][0] = 1

    radL[0][1] = 0
    midL[0][1] = 1

    radL[1][0] = 0.5*(1/X.a[1] - 1/X.b[1])
    midL[1][0] =  1/X.b[1] + radL[1][0]

    
    interval = -X[0]/(X[1]**2)
    radL[1][1] = 0.5*(interval.b - interval.a)
    midL[1][1] = interval.a + radL[1][1]
    return ip.Interval(midL, radL, midRadQ=True)

def F_x(x):
    B = ip.Interval([[3, 4], [1, 1]])
    F = ip.Interval([[0, 0], [0, 0]])

    F[0] = x[0] + x[1] - B[0]
    F[1] = x[0]/x[1] - B[1]
    return F

def inv_midA(A):
    midA = np.zeros((2, 2))
    for i in range (0,2):
        for j in range (0,2):
            midA[i][j] = A[i][j].mid
    return np.linalg.inv(midA)

def mid_X(X):
    midX = np.zeros(2)
    midX[0] = X[0].mid
    midX[1] = X[1].mid
    return midX

def matrix_norm(A):
    norm = 0
    sum = 0
    for i in range(0,2):
        sum = A[i][0].b + A[i][1].b
        if(sum > norm):
            norm = sum
    return norm

def b_norm(b):
    norm = 0
    num  = 0
    for i in range(0,2):
        num = b[i].b
        if num > norm:
            norm = num
    return norm


def dist(X, Y):
    d = np.zeros(2)
    for i in range(0,2):
        if(max(abs(X[i].b - Y[i].b), abs(X[i].a - Y[i].a)) < 0.0000001):
            return True
    return False

def Dist(X, Y):
    d = np.zeros(2)
    for i in range(0,2):
        d[i] = max(abs(X[i].b - Y[i].b), abs(X[i].a - Y[i].a))
    return d

def intersection(X, kr):
    midL = np.zeros(2)
    radL = np.zeros(2)
    a1 = kr.a[0] if (kr.a[0] >= X.a[0]) else X.a[0]
    b1 = kr.b[0] if (kr.b[0] <= X.b[0]) else X.b[0]
    a2 = kr.a[1] if (kr.a[1] >= X.a[1]) else X.a[1]
    b2 = kr.b[1] if (kr.b[1] <= X.b[1]) else X.b[1]

    radL[0] = (b1 - a1)/2
    midL[0] = a1 + radL[0]
    radL[1] = (b2 - a2)/2
    midL[1] = a2 + radL[1]

    return ip.Interval(midL, radL, midRadQ=True)


if __name__ == "__main__":
    midX = [1.75, 1.75]
    radX = [0.75, 0.75]
    I = [[1,0], [0,1]]
    X = ip.Interval(midX, radX, midRadQ=True)
    X_k = []
    i = 0
    print(X)
    while(True):
        print(i)
        print(X)
        i+=1
        X_k.append(X)
        A = create_L(X)
        x_av = mid_X(X)
        b = F_x(x_av)
        lambd = inv_midA(A)
        eta = matrix_norm(I - lambd @ A)
        teta = b_norm(abs(lambd @ b))/(1 - eta)
        # print(X,A, b,lambd,eta,teta, sep="\n")
        mid_kr = [0, 0]
        rad_kr = [teta, teta]
        X_kr = ip.Interval(mid_kr, rad_kr, midRadQ=True)
        kr = (lambd @ b + (I - lambd @ A) @ X_kr)
        new_X_kr = intersection(kr, X_kr)
        j = 0
        while(not dist(X_kr, new_X_kr)):
            j+=1
            X_kr = new_X_kr
            kr = (lambd @ b + (I - lambd @ A) @ X_kr)
            new_X_kr = intersection(kr, X_kr)
            # print(j, Dist(X_kr, new_X_kr) ,new_X_kr)
        X_kr = new_X_kr 
        N = x_av - X_kr
        print("N:", N)
        newX = intersection(X, N)
        print(Dist(X, newX))
        if(dist(X, newX)):
            X = newX
            break
        X = newX
    print(i, newX)

    plt.figure(1)
    x = np.arange(-10, 10, 0.01)
    y1_1 = 4 - x
    plt.plot(x, y1_1, color="blue")
    y1_2 = 3 - x
    plt.plot(x, y1_2, color="blue")
    plt.fill_between(y1=y1_1, y2=y1_2, x=x, color="blue",alpha=0.2,label=r"$x_1+x_2 = [3,4]$")
    y2 = x
    plt.plot(x, y2, color ="orange",label = r'$\frac{x_1}{x_2} = [1, 1]$')

    plt.grid()
    plt.xlim(-2, 6)
    plt.ylim(-2, 6)
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(x, y1_1, color="blue")
    plt.plot(x, y1_2, color="blue")
    plt.fill_between(y1=y1_1, y2=y1_2, x=x, color="blue",alpha=0.2,label=r"$x_1+x_2 = [3,4]$")
    plt.plot(x, y2, color ="orange",label = r'$\frac{x_1}{x_2} = [1, 1]$')
    one = abs(X_k[0][0].b - X_k[0][0].a)
    two = abs(X_k[0][1].b - X_k[0][1].a)
    mid1 = abs((X_k[0][0].a + X_k[0][0].b)/2)
    mid2 = abs((X_k[0][1].a + X_k[0][1].b)/2)
    print(one, two, mid1, mid2)
    iveRect = plt.Rectangle((X_k[0][0].a, X_k[0][1].a),one , two, edgecolor='black', facecolor='none', linewidth=2, label=r"$X^{(0)}$")
    plt.gca().add_patch(iveRect)
    plt.grid()
    plt.xlim(-1, 7)
    plt.ylim(-1, 7)
    plt.legend()
    plt.show()

    plt.figure(3)
    plt.plot(x, y1_1, color="blue")
    plt.plot(x, y1_2, color="blue")
    plt.fill_between(y1=y1_1, y2=y1_2, x=x, color="blue",alpha=0.2,label=r"$x_1+x_2 = [3,4]$")
    plt.plot(x, y2, color ="orange",label = r'$\frac{x_1}{x_2} = [1, 1]$')
    print(len(X_k))
    for i in range(0, len(X_k)):
        one = abs(X_k[i][0].b - X_k[i][0].a)
        two = abs(X_k[i][1].b - X_k[i][1].a)
        mid1 = abs((X_k[i][0].a + X_k[i][0].b)/2)
        mid2 = abs((X_k[i][1].a + X_k[i][1].b)/2)
        print(one, two, mid1, mid2)
        iveRect = plt.Rectangle((X_k[i][0].a, X_k[i][1].a), one , two, edgecolor='black', facecolor='none', linewidth=1.5)
        plt.gca().add_patch(iveRect)
    plt.grid()
    plt.xlim(0.5, 3)
    plt.ylim(0.5, 3)
    plt.legend()
    plt.show()

