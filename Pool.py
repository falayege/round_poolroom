import math 
import matplotlib.pyplot as plt
import random
import numpy as np

from Ball import Ball
from Hole import Hole

class Pool:
    def __init__(self, radius, friction, friction_rolling,friction_edge,rad,mu):
        """Initialisation of the board with its size and friction coefficient"""
        self.radius = radius
        self.friction = friction
        self.friction_rolling = friction_rolling
        self.friction_edge = friction_edge
        self.balls = []
        self.balls_in_motion = []
        self.holes = []
        self.white_ball = None
        self.hide_white_ball = False
        self.a_f= 9.81 * self.friction
        self.speed_stop = 0.2
        self.rho = 1.2
        self.C_d = 1.2
        self.mu = mu

        #for grid method
        self.grid_size = 2 * rad  # Assuming 'radius' is the radius of the largest ball
        self.grid = {} 


    def add_ball(self, ball: Ball) -> bool:
        """
        Add a ball to the board. 
        Returns True if the ball is added successfully, otherwise returns False.
        """
        # Check for overlap with existing balls
        for existing_ball in self.balls:
            dx = existing_ball.x - ball.x
            dy = existing_ball.y - ball.y
            distance = math.sqrt(dx**2 + dy**2)

            # If distance between centers is less than the sum of the radii, they overlap
            if distance < existing_ball.radius + ball.radius:
                return False
        if self.ball_in_hole(ball):
            return False
        self.balls.append(ball)
        if ball.speed()>0:
            self.balls_in_motion.append(ball)
        # Update the grid with the new ball
        self.add_ball_to_grid(ball)
        return True
    
    def add_ball_to_grid(self, ball: Ball):
        """
        Add ball to grid
        """
        grid_row = int(ball.y / self.grid_size)
        grid_col = int(ball.x / self.grid_size)
        grid_key = (grid_row, grid_col)
        if grid_key not in self.grid:
            self.grid[grid_key] = []
        self.grid[grid_key].append(ball)
        ball.grid_cell = grid_key

    def remove_ball_from_grid(self, ball: Ball):
        if ball.grid_cell in self.grid and ball in self.grid[ball.grid_cell]:
            self.grid[ball.grid_cell].remove(ball)

    def add_rand_ball(self,radius,magnus_cst,color=None)->bool:
        """
        Add a ball generated randomly
        If the ball can be added to the board, it's added and the method return True
        """
        angle = random.uniform(-math.pi, math.pi)  
        max_distance = self.radius - radius  
        distance = random.uniform(0, max_distance) 
        x = distance * math.cos(angle) 
        y = distance * math.sin(angle) 
        ball = Ball(x, y, radius,magnus_cst,0,0,color) 
        return(self.add_ball(ball))
    
    def add_white_black_balls(self,radius,magnus_cst):
        """
        Add the two balls of pool the black and the white
        """
        while not(self.add_rand_ball(radius,magnus_cst,'white')):
            pass
        self.white_ball = self.balls[0]
        while not(self.add_rand_ball(radius,magnus_cst,'black')):
            pass
        """
        self.add_ball(Ball(0.5,0.5,radius,0,0,magnus_cst,'black'))
        self.add_ball(Ball(0.2,0.2,radius,5,4,magnus_cst,'white'))
        """

    def remove_ball(self, ball: Ball):
        "Remove a ball"
        ball.vx = 0
        ball.vy = 0
        if ball in self.balls:
            if ball is self.white_ball:  # If the removed ball is the white ball
                self.hide_white_ball = True
                # Display a banner for the white ball
                plt.gcf().text(0.5, 0.5, 'White ball potted!', 
                            horizontalalignment='center', verticalalignment='center', 
                            fontsize=14, color='white', bbox=dict(facecolor='red', alpha=0.8))
                plt.pause(1)  # Display the message for 1 second
                plt.gcf().texts.clear()  # Clear the text
                plt.draw()  # Redraw the current figure to update the display
            self.remove_ball_from_grid(ball)
            self.balls.remove(ball)
        if ball in self.balls_in_motion:
            self.balls_in_motion.remove(ball)

    def add_hole(self,hole: Hole):
        """Add a hole to the board"""
        self.holes.append(hole)

    def add_usual_holes(self,rayon):
        """Add holes"""
        self.holes.append(Hole(math.sqrt(2)*self.radius/2, math.sqrt(2)*self.radius/2, rayon))
        self.holes.append(Hole(math.sqrt(2)*self.radius/2, -math.sqrt(2)*self.radius/2, rayon))
        self.holes.append(Hole(-math.sqrt(2)*self.radius/2, -math.sqrt(2)*self.radius/2, rayon))
        self.holes.append(Hole(-math.sqrt(2)*self.radius/2, math.sqrt(2)*self.radius/2, rayon))


    def ball_in_hole(self,ball:Ball)->bool:
        "chek if a ball is in a hole"
        for hole in self.holes:
            dx = hole.x - ball.x
            dy = hole.y - ball.y
            distance = math.sqrt(dx**2 + dy**2)
            if distance < ball.radius/4+hole.radius: #consider that it falls into the hole
                return True
        return False
    
    def upgrade_loc_hole(self,dt):
        for hole in self.holes:
            angle = math.atan2(hole.y, hole.x)
            angle = (angle + self.mu * dt) % (2 * math.pi)
            hole.x = math.cos(angle) * self.radius
            hole.y = math.sin(angle) * self.radius

    def all_balls_stopped(self):
        """Check if all balls have come to a stop."""
        for ball in self.balls:
            if ball.speed() > 0.2:  # Using the threshold you've defined for "stopped" 
                return False
        return True
    
    def check_collisions(self, ball: Ball):
        """Check collisions for a given ball with balls in adjacent cells and in self cell"""
        grid_row, grid_col = ball.grid_cell
        for row in range(grid_row - 1, grid_row + 2):
            for col in range(grid_col - 1, grid_col + 2):
                neighbor_key = (row, col)
                if neighbor_key in self.grid:
                    for other_ball in self.grid[neighbor_key]:
                        if ball != other_ball:
                            ball.collide(other_ball)

    def get_adjacent_cells_to_motion(self):
        adjacent_cells = set()
        for ball in self.balls_in_motion:
            grid_row, grid_col = ball.grid_cell
            for row in range(grid_row - 1, grid_row + 2):
                for col in range(grid_col - 1, grid_col + 2):
                    adjacent_cells.add((row, col))
        return adjacent_cells
    
    def get_distance(self, x1, y1, x2=0, y2=0):
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    def euler_explicit(self,ball, dt):
        net_force = self.calculate_forces(ball,ball.vx,ball.vy)
        ball.vx += net_force[0] * dt
        ball.vy += net_force[1] * dt
        ball.move(dt)

    def runge_kutta_4(self, ball, dt):
        # K1
        k1_v = self.calculate_forces(ball, ball.vx, ball.vy)
        k1_p = np.array([ball.vx, ball.vy])

        # K2
        k2_v = self.calculate_forces(ball, ball.vx + dt/2 * k1_v[0], ball.vy + dt/2 * k1_v[1])
        k2_p = np.array([ball.vx, ball.vy]) + dt/2 * k1_v

        # K3
        k3_v = self.calculate_forces(ball, ball.vx + dt/2 * k2_v[0], ball.vy + dt/2 * k2_v[1])
        k3_p = np.array([ball.vx, ball.vy]) + dt/2 * k2_v

        # K4
        k4_v = self.calculate_forces(ball, ball.vx + dt * k3_v[0], ball.vy + dt * k3_v[1])
        k4_p = np.array([ball.vx, ball.vy]) + dt * k3_v

        # Update the velocity
        ball.vx += dt/6 * (k1_v[0] + 2*k2_v[0] + 2*k3_v[0] + k4_v[0])
        ball.vy += dt/6 * (k1_v[1] + 2*k2_v[1] + 2*k3_v[1] + k4_v[1])

        # Update the position
        ball.x += dt/6 * (k1_p[0] + 2*k2_p[0] + 2*k3_p[0] + k4_p[0])
        ball.y += dt/6 * (k1_p[1] + 2*k2_p[1] + 2*k3_p[1] + k4_p[1])

    def rk2_step(self, ball, dt):
        # Étape 1: Calculer k1 (force initiale)
        k1_v = self.calculate_forces(ball, ball.vx, ball.vy) * dt
        k1_p = np.array([ball.vx, ball.vy]) * dt

        # Étape 2: Calculer k2 en utilisant k1 pour prédire la position et la vitesse au milieu de l'intervalle
        k2_v = self.calculate_forces(ball, ball.vx + k1_v[0]/2, ball.vy + k1_v[1]/2) * dt
        k2_p = (np.array([ball.vx, ball.vy]) + k1_v/2) * dt

        # Mise à jour finale en utilisant k2
        ball.vx += k2_v[0]
        ball.vy += k2_v[1]
        ball.x += k2_p[0]
        ball.y += k2_p[1]

    def update_implicit(self, ball, dt):
        # Estimation de la position et de la vitesse futures
        future_vx = ball.vx + self.calculate_forces(ball, ball.vx, ball.vy)[0] * dt
        future_vy = ball.vy + self.calculate_forces(ball, ball.vx, ball.vy)[1] * dt
        future_x = ball.x + future_vx * dt
        future_y = ball.y + future_vy * dt

        # forces calcuated on the based speed/positions estimates
        net_force = self.calculate_forces(ball, future_vx, future_vy)

        # update position and speed
        ball.vx += net_force[0] * dt
        ball.vy += net_force[1] * dt
        ball.x = future_x
        ball.y = future_y


    def calculate_forces(self, ball, vx, vy):
        """
        Calculate all the forces applied on a ball at a given moment
        """
        if ball.is_sliding():   
            # Coulomb friction
            if np.linalg.norm([vx, vy]) > 0:
                friction_direction = -np.array([vx, vy]) / np.linalg.norm([vx, vy])
                friction_force = self.a_f * friction_direction
            else:
                friction_force = np.array([0, 0])
        else:
            # Rolling friction
            if np.linalg.norm([vx, vy]) > 0:
                rolling_friction_magnitude = self.friction_rolling*ball.rolling_speed()**2
                friction_direction = -np.array([vx, vy]) / np.linalg.norm([vx, vy])
                friction_force =  rolling_friction_magnitude* friction_direction
            else:
                friction_force = np.array([0, 0])

        # Air resistance
        v = np.linalg.norm([vx, vy])
        if v != 0:
            drag_force_magnitude = 0.5 * self.C_d * self.rho * math.pi * (ball.radius ** 2) * (v ** 2)/ball.mass
            drag_force = drag_force_magnitude * np.array([vx, vy]) / v
        else:
            drag_force = np.array([0, 0])

        # Coriolis effect
        coriolis_force = 2 * np.cross([vx, vy, 0], [0, 0, self.mu])

        #Magnus effect
        if ball.omega != 0:
            magnus_effect = ball.magnus_cst * ball.omega
            magnus_force = np.array([-magnus_effect * vy, magnus_effect * vx])
        else:
            magnus_force = np.array([0, 0])

        return (friction_force + drag_force + coriolis_force[:2]+magnus_force)

    def transform_to_inertial_reference(self, dt):
        # For each ball, transform its position and velocity from the table frame to the inertial frame
        for ball in self.balls:                
            angular_displacement = self.mu * dt
            distance_to_center = math.sqrt(ball.x**2 + ball.y**2)
            angle_to_center = math.atan2(ball.y, ball.x)
            new_angle_to_center = angle_to_center + angular_displacement

            # Update ball position in the inertial frame
            ball.x = distance_to_center * math.cos(new_angle_to_center)
            ball.y = distance_to_center * math.sin(new_angle_to_center)


    def step(self, dt):
        self.upgrade_loc_hole(dt)
        """Update of the position of the balls  and thier speed during a dt interval"""
        update_balls_in_motion = set()  # Using a set to ensure uniqueness
        for ball in self.balls_in_motion:
            self.remove_ball_from_grid(ball)
            if self.ball_in_hole(ball):
                if isinstance(ball.color,str):
                    if ball.color == 'black':
                        plt.text(self.radius / 2, self.radius / 2, 'Game Over! Black ball fell into a hole.', 
                                horizontalalignment='center', verticalalignment='center', 
                                fontsize=14, bbox=dict(facecolor='red', alpha=0.5))
                        plt.pause(2)  # Pause to display the message before closing
                        plt.close()
                        print("Game Lost! Black ball fell into a hole.")
                        exit()
                    elif ball.color == 'white':
                        angle = random.uniform(0, 2 * math.pi)  # Random angle in radians
                        max_distance = self.radius - ball.radius  # Maximum distance from the center
                        distance = random.uniform(0, max_distance)  # Random distance from the center
                        ball.x = distance * math.cos(angle)  # X coordinate
                        ball.y = distance * math.sin(angle)  # Y 
                        ball.vx = 0  # Reset velocity
                        ball.vy = 0
                        ball.omega = 0
                        self.remove_ball(ball)  # Remove the white ball from the game
                    continue
                else:
                    self.remove_ball(ball)
                    continue
            self.euler_explicit(ball,dt)
            self.add_ball_to_grid(ball)
            self.check_collisions(ball)                    


            distance_to_center = self.get_distance(ball.x, ball.y)

            # Check if the ball is colliding with the edge
            if distance_to_center + ball.radius > self.radius:
                # Normalize the direction vector from the ball to the table's center
                dx = (ball.x - 0) / distance_to_center
                dy = (ball.y - 0) / distance_to_center
                ball.x*=(distance_to_center-ball.radius/6)/distance_to_center
                ball.y*=(distance_to_center-ball.radius/6)/distance_to_center
                ball.rebound(dx, dy, self.friction_edge)

        # Calculate the rebound direction (reflecting the direction vector about the normal)

            # If ball is still moving, add to the list of balls in motion
            if ball.speed() > self.speed_stop:  # Consider a threshold for "stopped"
                update_balls_in_motion.add(ball)

        if self.mu !=0:
            for ball in self.balls:
                if ball not in self.balls_in_motion:
                    self.remove_ball_from_grid(ball)
                    self.add_ball_to_grid(ball)
            
        # Adding the balls in motion to the list
        for ball in self.balls:
            if ball.speed() > self.speed_stop: 
                update_balls_in_motion.add(ball)
            else:
                ball.omega = 0
        
            
        # Update the list of balls in motion
        self.balls_in_motion = list(update_balls_in_motion)

        if self.all_balls_stopped():
        #If the white ball should be added back into the game - once fallen in a hole
            if self.hide_white_ball:
                self.balls.insert(0,self.white_ball) #insert the white ball at its position
                self.hide_white_ball = False
                self.balls_in_motion.append(self.white_ball)
                self.add_ball_to_grid(self.white_ball)

        if self.mu != 0:
            self.transform_to_inertial_reference(dt)


