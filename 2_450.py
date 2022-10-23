import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

# u(0, t) = F_0
# u(l, t) = 0

#Task parameters
a = 1
r_0 = 1
C = 0.001
alpha = 1

#Параметры сетки
A = 0
B = r_0
N = 200

r = np.linspace(A, B, N)

# Выводные параметры задачи
h  = (B - A)/(N - 1)
dt = h ** 2/ 2 / a**2

# Время
time     = 0.0
end_time = 1000.0

# Инициализация u1, u2
u1 = np.zeros(N)


def solution_step():
    global time, u1
    
    for _ in range(100):
        # Основная часть
        u2 = np.zeros(N)
        for i in range(1, N - 1):
            u2[i] = u1[i] +  dt*a**2 * (2/r[i]*(u1[i] - u1[i - 1])/h + (u1[i+1] - 2*u1[i] + u1[i - 1])/h**2)
        # Граничные условия
        u2[0]  = (- u2[2] + 4 * u2[1]) / 3
        u2[N - 1]  =  u2[N - 2] + h*alpha*time**2

        u1 = u2
        time += dt


fig = plt.figure(1, dpi=200)
ax = plt.axes(xlim=(A, B) , ylim=(-0.5, 0.5))
line, = plt.plot(r, u1)
ax.set_xlabel("x")
ax.set_ylabel("u")
title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    title.set_text("t = %.3f" % time)
    return line, title, 


def animate(frame):
    global u1, u2, time
    if time > end_time:
        raise 'Enough'

    solution_step()
    line.set_data(r, u1)
    #ax.set_ylim(ymax=np.max(u1), ymin=np.min(u1))
    title.set_text("t = %.3f" % time)
    return line, title, 

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=10, blit=True)

plt.show()
