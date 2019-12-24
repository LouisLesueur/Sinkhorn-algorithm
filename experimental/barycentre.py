import numpy as np
import time
import matplotlib.pyplot as plt
from cv2 import cv2



class SimplexFromImage:
    def __init__(self, image):
        self.img_size = len(image)
        self.lenght = len(image)**2

        self.simplex = [image[i, j] for i in range(
            self.img_size) for j in range(self.img_size)]
        self.constant = sum(self.simplex)
        self.values = (1/self.constant)*np.array(self.simplex)


class ImageFromSimplex:
    def __init__(self, simplex, img_size):
        self.img_size = img_size
        self.pixels = np.zeros((self.img_size, self.img_size))
        for i in range(self.img_size):
            for j in range(self.img_size):
                self.pixels[i, j] = simplex[i*self.img_size+j]
    def show(self):
        plt.imshow(self.pixels, cmap="gray")
        plt.show()




def barycenter(p1, p2, lamb, eps, n_iter):
    m = p1.img_size
    m2 = m**2
    ksi = np.zeros((m2, m2))
    for i in range(m2):
        for j in range(m2):
            x_i = (i+1/2)/m2
            x_j = (j+1/2)/m2
            dis = (x_i-x_j)**2
            ksi[i, j] = np.exp(-dis/eps)
    u1 = u2 = np.ones(m2)
    for i in range(n_iter):
        # Computing vk(n+1)
        v1 = np.divide(p1.values, np.dot(ksi.T, u1))
        v2 = np.divide(p2.values, np.dot(ksi.T, u2))
        # Computing p(n+1)
        first_term = np.power(np.dot(ksi.T, v1), lamb)
        second_term = np.power(np.dot(ksi.T, v2), 1-lamb)
        p = np.multiply(first_term, second_term)
        # Computing uk(n+1)
        u1 = np.divide(p, np.dot(ksi.T, v1))
        u2 = np.divide(p, np.dot(ksi.T, v2))
    # Computing projections matrix
    gamma1 = np.dot(np.diag(u1), np.dot(ksi, np.diag(v1)))
    gamma2 = np.dot(np.diag(u2), np.dot(ksi, np.diag(v2)))
    q1 = p1.constant * np.dot(gamma1.T, np.ones(m2))
    q2 = p2.constant * np.dot(gamma2.T, np.ones(m2))
    bar1 = (lamb*p1.constant + (1-lamb)*p2.constant)*np.dot(gamma1, np.ones(m2))
    bar2 = (lamb*p1.constant + (1-lamb)*p2.constant)*np.dot(gamma2, np.ones(m2))
    print(max(bar1))
    print(max(bar2))
    return q1, q2, bar1, bar2, p


def barycenter2(p1, p2, lamb, eps, n_iter):
    lamb1 = lamb
    lamb2 = 1 - lamb
    m = p1.img_size**2
    K = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            x_i = (i+1/2)/m
            x_j = (j+1/2)/m
            dis = (x_i-x_j)**2
            K[i, j] = np.exp(-dis/eps)
    a1 = a2 = np.ones(m)
    for l in range(n_iter):
        # computing p(l)
        p = np.multiply(np.power(np.dot(K.T,a1), lamb1), np.power(np.dot(K.T, a2), lamb2))
        # computing bs(l+1)
        b1 = np.divide(p, np.dot(K.T, a1))
        b2 = np.divide(p, np.dot(K.T, a2))
        # computing as(l+1)
        a1 = np.divide(p1.values, np.dot(K, b1))
        a2 = np.divide(p2.values, np.dot(K, b2))
    return (lamb1*p1.constant + lamb2*p2.constant)*p

def barycenter2_curves(y1, y2, lamb, eps, n_iter):
    lamb1 = lamb
    lamb2 = 1 - lamb
    m = len(y1)
    K = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            x_i = (i+1/2)/m
            x_j = (j+1/2)/m
            dis = (x_i-x_j)**2
            K[i, j] = np.exp(-dis/eps)
    a1 = a2 = np.ones(m)
    for l in range(n_iter):
        # computing p(l)
        p = np.multiply(np.power(np.dot(K.T,a1), lamb1), np.power(np.dot(K.T, a2), lamb2))
        # computing bs(l+1)
        b1 = np.divide(p, np.dot(K.T, a1))
        b2 = np.divide(p, np.dot(K.T, a2))
        # computing as(l+1)
        a1 = np.divide(y1, np.dot(K, b1))
        a2 = np.divide(y2, np.dot(K, b2))
    return p

T0 = time.time()

SQUARE = cv2.imread('./SQUARE_32.png', 0)
CIRCLE = cv2.imread('./CIRCLE_32.png', 0)

p1 = SimplexFromImage(SQUARE)
p2 = SimplexFromImage(CIRCLE)

lamb = 0.8
eps = 2/(p1.lenght)
n_iter = 100
print("temps écoulé = " + str(time.time()-T0))


X = np.linspace(0,1,100)
Y1 = [np.exp(x) for x in X]
Y1 = Y1/sum(Y1)
Y2 = [np.cos(2*np.pi*x)+1 for x in X]
Y2 = Y2/sum(Y2)


plt.figure(1)
plt.plot(X, Y1, 'b')
plt.plot(X, Y2, 'r')
bar = barycenter2_curves(Y1, Y2, lamb, eps, n_iter)
plt.plot(X, bar, 'g')
plt.show()
