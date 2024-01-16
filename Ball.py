import math
import random
import matplotlib.colors as mcolors
import numpy as np

class Ball:
    """
    The class Ball represents a ball on a poolroom
    A ball in defined by its position and speed in 2D
    """
    def __init__(self, x, y, radius, magnus_cst, vx=0, vy=0,color=None,omega=0):
        """Initialization of a ball with its position and speed."""
        self.x = x
        self.y = y
        self.radius = radius
        self.w_max = 6/self.radius #6m/s
        self.vx = vx
        self.vy = vy
        self.mass = 4/3*math.pi*radius**3*2000
        self.magnus_cst = magnus_cst
        self.omega = omega  # Angular velocity
        self.grid_cel = None
        self.marker_angle = 0 #initial angle of the marker
        if color is None:
            self.color = mcolors.hsv_to_rgb([np.random.rand(), 1, 1])
        else:
            self.color = color

        
    def move(self, dt):
        """Move the ball during an interval dt"""
        self.x += self.vx * dt
        self.y += self.vy * dt
        
    def magnus(self,dt):
        """
        Handle the magnus effect on the ball
        """
        magnus_effect = self.magnus_cst * self.omega
        self.vx += magnus_effect * self.vy * dt
        self.vy -= magnus_effect * self.vx * dt


    def speed(self):
        return math.sqrt(self.vx**2 + self.vy**2)
    
    def rolling_speed(self):
        return self.speed()/self.radius
    
    def is_sliding(self):
        return self.rolling_speed()>self.w_max

    def rebound(self, normal_x, normal_y,friction_edge):
        """Rebound of the ball"""
        dot = self.vx * normal_x + self.vy * normal_y
        self.vx -= 2 * dot * normal_x
        self.vx *= (1-friction_edge)
        self.vy -= 2 * dot * normal_y
        self.vy *= (1-friction_edge)
        self.omega *=-(1-friction_edge)
        
    
    def collide(self, other):
        """
        Handle collision with another ball other
        Update the speed of the two balls
        Collision between the balls is elastic no loss of cinetic energy
        """
        dx = self.x - other.x
        dy = self.y - other.y
        distance_pow_2 = dx**2 + dy**2

        if (self.radius + other.radius)**2 < distance_pow_2: # check if collision
            return

        dvx = self.vx - other.vx
        dvy = self.vy - other.vy

        dot = dvx * dx + dvy * dy

        if dot  > 0:  # for the balls not to be stuck on
            return
        
        impact = dot/distance_pow_2
    
        self.vx -= impact * dx
        self.vy -= impact * dy
        other.vx += impact * dx
        other.vy +=  impact * dy

        # Normalized direction vector from other to self
        nx = dx / distance_pow_2**.5
        ny = dy / distance_pow_2**.5

        # Collision point coordinates
        collision_x = other.x + nx * other.radius
        collision_y = other.y + ny * other.radius

        # Calculate the angle of impact
        impact_dx = collision_x - self.x
        impact_dy = collision_y - self.y
        impact_distance = math.sqrt(impact_dx**2 + impact_dy**2)

        if impact_distance < 1e-6:  # Avoid division by zero
            impact_angle = 0
        else:
            cos_impact_angle = (impact_dx * nx + impact_dy * ny) / impact_distance
            cos_impact_angle = max(-1, min(1, cos_impact_angle))  # Clamping
            impact_angle = math.acos(cos_impact_angle)
        # Determine if it's a center hit based on the impact angle
        threshold_angle = math.radians(10)
        center_collision_self = impact_angle < threshold_angle
        if not center_collision_self:
            self.omega += dot / self.radius

        other_nx = -nx 
        other_ny = -ny

        # Calculate the angle of impact for the other ball
        other_impact_dx = collision_x - other.x
        other_impact_dy = collision_y - other.y
        other_impact_distance = math.sqrt(other_impact_dx**2 + other_impact_dy**2)

        if other_impact_distance < 1e-6:  # Avoid division by zero
            other_impact_angle = 0
        else:
            # Calculate cosine of the angle and clamp it between -1 and 1
            cos_other_impact_angle = (other_impact_dx * other_nx + other_impact_dy * other_ny) / other_impact_distance
            cos_other_impact_angle = max(-1, min(1, cos_other_impact_angle)) 

            other_impact_angle = math.acos(cos_other_impact_angle)
                
        center_collision_other = other_impact_angle < threshold_angle

        if not center_collision_other:
            other.omega -= dot / other.radius