import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

# u(0, t) = F_0
# u(l, t) = 0

#Task parameters
a = 1
l = 10
T = 1 # Сила натяжения струны
F_0 = 0.01 # Начальная сила, действующая на струну до t=0

C = 0.95

#Параметры сетки
A = 0.0
B = l
N = 301

x = np.linspace(A, B, N)

# Выводные параметры задачи
h  = (B - A)/(N - 1)
dt = C * h / a

# Время
time     = 0.0
end_time = 1000.0

# Инициализация u1, u2
u1 = np.zeros(N) 
u2 = F_0/T*((x - l)/l) # Начальное условие

def solution_step():
    global time, u1, u2
    
    for _ in range(10):
        # Основная часть
        u3 = np.zeros(N)
        for i in range(1, N - 1):
            u3[i] = 2 * u2[i] - u1[i] + (a * dt / h) ** 2 * (u2[i + 1] - 2 * u2[i] + u2[i - 1])
        # Граничные условия
        u3[0]  = (- u3[2] + 4 * u3[1]) / 3 # Взято из numeric'a
        u3[N - 1] = 0

        u1, u2 = u2, u3

        time += dt

fig = plt.figure(1, dpi=200)
ax = plt.axes(xlim=(A, B), ylim=(-3, 3))
line, = plt.plot(x, u1)
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
    line.set_data(x, u1)
    #ax.set_ylim(ymax=np.max(u1), ymin=np.min(u1))
    title.set_text("t = %.3f" % time)
    return line, title, 

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=10, blit=True)

plt.show()
