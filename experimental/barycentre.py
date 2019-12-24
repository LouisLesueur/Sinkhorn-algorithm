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

T0 = time.time()

SQUARE = cv2.imread('./SQUARE_32.png', 0)
CIRCLE = cv2.imread('./CIRCLE_32.png', 0)

p1 = SimplexFromImage(SQUARE)
p2 = SimplexFromImage(CIRCLE)

lamb = 1
eps = 2/p1.lenght
n_iter = 100
q1, q2, bar1, bar2, p = barycenter(p1, p2, lamb, eps, n_iter)
print("temps écoulé = " + str(time.time()-T0))
q1 = ImageFromSimplex(q1, p1.img_size)
q1.show()
q2 = ImageFromSimplex(q2, p2.img_size)
q2.show()
bar1 = ImageFromSimplex(bar1, p1.img_size)
bar1.show()
bar2 = ImageFromSimplex(bar2, p1.img_size)
bar2.show()
print("aaa")
bar3 = ImageFromSimplex(p, p1.img_size)
bar3.show()
