# Freekick-simulation
The code file numerically simulates the tracjectory of a football during a free kick, including gravity, air drag, magnus effect depending on wind conditions.
The blue line in the animation shows how the football would move after the ball is kicked in a football field.
To carry out the simulation, you can adjust the parameters in the physical constraint and initial condition session. Details are given for coressponding variables.  
The point where the ball is kicked from the centre could be adjusted as it would induce a spin and contribute with the magnus effect:
```
kick_pt = np.array([-0.1, 0.01, 0])
```
You can define the wind speed as a vector field, which could be dependant in position or time
```
def wind(pos, t):
    x, y, z = pos
    wind_x = 5+0.5 * np.sin(0.2 * t + 0.1 * y)
    wind_y = -3+ 0.2 * np.cos(0.3 * t + 0.1 * x)  # Push rightward
    wind_z = 0

    return np.array([wind_x, wind_y, wind_z])
```
The posistion where the free kick is initiated can also be defined as below
```
pos = [np.array([15, 5, 0.2])]
```
After setting all the parameters you can whether or not to make the 3D graph in equal scale by calling the "equal_scale" function:
```
equal_scale(ax)
```
Finally, you can enter the following commond to the terminal and run the simulation:
```
python3 freekick_path.py
```

# How the physics and calculation works
Before explaining how the simulation is built, a few assumptions are made in the simulation:
- Ball is a perfect sphere
- Wind conditions are in laminar flow
- Air density is constant
- Gravity variation is neglected
- The part where the ball is kicked is modelled as a point
- Spin axis remains constant during the flight
- No contact or deformation during the flight
Whatsmore, the simulation terminates when the ball reaches the ground or goes out of bounds to minimise randomness during the simulation such as reaction of goal keeper or players.
#
This simulation includes all possible forces during a free kick such as gravity, air drag and magnus force. The magnus effect comes from the varible kick_pt, which is an offset vector from the ball’s center to the point of contact. v0_vec is the initial velocity vector (direction you’re kicking in). The cross product gives you the spin axis (omega), which shows the direction and magnitude of rotation. This rotation causes a difference in air flow on both sides. Benoulli's principal shows that there would be a pressure difference due to the varying flow speed, inducing a magnus force. The other forces are determined by relevant equations. In each time step, these forces are calculated depending on the current posisiton of the ball and the net force is found by newton's second law. The acceleration value is yielded and RK4 integration method is carried out to determine the position and velocity for hte next time step. The loop breaks when the ball goes out of bounds or touches the ground. Finally, the data points are saved and ploted with animation to show the simulation results. 
