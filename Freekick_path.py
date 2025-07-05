
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
# Physical constraints:
m = 0.43  # mass (kg)
r = 0.11  # radius (m)
A = np.pi * r**2 # cross-sectional area
rho = 1.2  # air density
Cd = 0.12 # drag coefficient
Cl = 2 # lift coefficient (for Magnus effect)
g = 9.81 #gravity

# Goal and field dimensions
goal_pos = 40
goal_width = 7.32
goal_height = 2.44
field_length = 40
field_width = 30
left_bot = [goal_pos, -goal_width / 2, 0]
right_bot = [goal_pos, goal_width / 2, 0]
left_top = [goal_pos, -goal_width / 2, goal_height]
right_top = [goal_pos, goal_width / 2, goal_height]

# Initial condition
dt = 0.01
t_max = 4
n = int(t_max / dt) #time steps
v0 = 28  # inital speed
theta = np.radians(10)  # vertical angle
phi = np.radians(5)    # horizontal angle

# Initial velocity vector components
vx0 = v0 * np.cos(theta) * np.cos(phi) 
vy0 = v0 * np.cos(theta) * np.sin(phi) 
vz0 = v0 * np.sin(theta)
v0_vec = np.array([vx0, vy0, vz0]) #velocity vector

#Rotation 
kick_pt = np.array([-0.1, 0.01, 0])  # spin direction source 
omega = np.cross(kick_pt, v0_vec)


# Restrictions to the initial kick point 
if np.linalg.norm(kick_pt) > 0.8 * r: #Arbitary restriction
    raise ValueError("Kick point too close to edge") # unrealistic for actual kick.


#Wind speed
def wind(pos, t):
    x, y, z = pos

    wind_x = 5+0.5 * np.sin(0.2 * t + 0.1 * y)
    wind_y = -3+ 0.2 * np.cos(0.3 * t + 0.1 * x)  # Push rightward
    wind_z = 0

    return np.array([wind_x, wind_y, wind_z])

# Magnus force function
def magnus_force(v, omega, v_wind):
    v_rel = v - v_wind
    v_mag = np.linalg.norm(v_rel)
    if v_mag == 0:
        return np.zeros(3)
    return 0.5 * rho * Cl * A * v_mag * np.cross(omega, v_rel / v_mag)# just the equation

# Drag force function (Due to air resistance)
def drag_force(v_rel):
    v_mag = np.linalg.norm(v_rel)
    if v_mag == 0:
        return np.zeros(3)
    return -0.5 * rho * Cd * A * v_mag * v_rel  #just the equation

# RK4 integration
pos = [np.array([15, 5, 0.2])] # initial position (x, y, z) for the free kick
vel = [v0_vec] #initial velocity vector
times = [0]

for _ in range(n):
    t = times[-1]
    x, v = pos[-1], vel[-1]

    def acceleration(v_local):
        v_wind = wind(x, t)
        f_drag = drag_force(v_local - v_wind)
        f_magnus = magnus_force(v_local, omega, v_wind)
        f_gravity = np.array([0, 0, -m * g])
        return (f_drag + f_magnus + f_gravity) / m

    k1v = acceleration(v) * dt
    k1x = v * dt

    k2v = acceleration(v + 0.5 * k1v) * dt
    k2x = (v + 0.5 * k1v) * dt

    k3v = acceleration(v + 0.5 * k2v) * dt
    k3x = (v + 0.5 * k2v) * dt

    k4v = acceleration(v + k3v) * dt
    k4x = (v + k3v) * dt

    new_v = v + (k1v + 2 * k2v + 2 * k3v + k4v) / 6
    new_x = x + (k1x + 2 * k2x + 2 * k3x + k4x) / 6

    if (
        new_x[2] < 0 or
        new_x[0] > field_length or
        abs(new_x[1]) > field_width / 2
    ):
        break

    pos.append(new_x)
    vel.append(new_v)
    times.append(t + dt)

pos = np.array(pos)

# Animation

# Draw goal
def draw_goal_and_field(ax):
    # Goalpost
    ax.plot([left_bot[0], right_bot[0]], [left_bot[1], right_bot[1]], [0, 0], 'r-')
    ax.plot([left_bot[0], left_top[0]], [left_bot[1], left_top[1]], [left_bot[2], left_top[2]], 'r-')
    ax.plot([right_bot[0], right_top[0]], [right_bot[1], right_top[1]], [right_bot[2], right_top[2]], 'r-')
    ax.plot([left_top[0], right_top[0]], [left_top[1], right_top[1]], [left_top[2], right_top[2]], 'r-')

    # Court lines
    ax.plot([0, 0], [-field_width / 2, field_width / 2], [0, 0], 'gray')
    ax.plot([0, field_length], [-field_width / 2, -field_width / 2], [0, 0], 'gray')
    ax.plot([0, field_length], [field_width / 2, field_width / 2], [0, 0], 'gray')
    ax.plot([field_length, field_length], [-field_width / 2, field_width / 2], [0, 0], 'red')

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([0, field_length])
ax.set_ylim([-field_width / 2, field_width / 2])
ax.set_zlim([0, max(np.max(pos[:, 2]), goal_height) + 1])
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Football Free Kick Simulation")

line, = ax.plot([], [], [], lw=2, color='blue' ,label ="Ball Path")
point, = ax.plot([], [], [], 'ro', label="Ball")

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    draw_goal_and_field(ax)
    ax.legend()
    return line, point

def animate(i):
    line.set_data(pos[:i, 0], pos[:i, 1])
    line.set_3d_properties(pos[:i, 2])
    point.set_data(pos[i, 0], pos[i, 1])
    point.set_3d_properties(pos[i, 2])
    return line, point

def equal_scale(ax, zmin=0):
    '''Set 3D plot axes to equal scale.'''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = np.ptp(x_limits)
    y_range = np.ptp(y_limits)
    z_range = np.ptp(z_limits)
    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_max = z_limits[1]
    z_middle = (z_max + zmin) / 2

    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([zmin, zmin + max_range])

ani = animation.FuncAnimation(fig, animate, frames=len(pos), init_func=init,
                              interval=30, blit=True)
equal_scale(ax) # only call if you want to set a actual scale
plt.tight_layout()
plt.show()