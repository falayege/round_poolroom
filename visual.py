import math
import random
import matplotlib.pyplot as plt
from Ball import Ball
from Pool import Pool
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as patches
from matplotlib.widgets import Slider, Button


def visualize_pool(pool,nb_balls,dt):
    """
    Allow the visual representation of the pool
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(-pool.radius, pool.radius)
    ax.set_ylim(-pool.radius, pool.radius)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    # Drawing the circular pool table
    table_circle = patches.Circle((0, 0), pool.radius, linewidth=2, edgecolor='brown', facecolor='darkgreen')
    ax.add_patch(table_circle)
    
    # Drawing holes
    hole_circles = {}
    for hole in pool.holes:
        hole_circle = patches.Circle((hole.x, hole.y), hole.radius, fc='black')
        plt.gca().add_patch(hole_circle)
        hole_circles[hole] = hole_circle


    # Initialize ball_to_circle dictionary
    ball_to_circle = {ball: plt.Circle((ball.x, ball.y), ball.radius, fc=ball.color) for ball in pool.balls}
    for circle in ball_to_circle.values():
        ax.add_patch(circle)

    axstop = plt.axes([0.77, 0.01, 0.1, 0.07])
    btn_stop = Button(axstop, 'Stop Game')

    def stop_game(event):
        print("Game Stopped")
        exit()
            
    btn_stop.on_clicked(stop_game)


    def update_velocity(val):
        if pool.all_balls_stopped():
            pool.balls[0].vx = vx_slider.val
            pool.balls[0].vy = vy_slider.val

    def update_rotation_speed(val):
        pool.mu = rotation_speed_slider.val


    axcolor = 'lightgoldenrodyellow'
    axvy = plt.axes([0.2, 0.04, 0.5, 0.03], facecolor=axcolor)
    axvx = plt.axes([0.2, 0.08, 0.5, 0.03], facecolor=axcolor)
    vx_slider = Slider(axvx, " White ball's speed (Vx)", -50.0, 50.0, valinit=pool.balls[0].vx)
    vy_slider = Slider(axvy, " White ball's speed (Vy)", -50.0, 50.0, valinit=pool.balls[0].vy)
    vx_slider.on_changed(update_velocity)
    vy_slider.on_changed(update_velocity)
    axcolor = 'lightgoldenrodyellow'
    axrotation = plt.axes([0.2, 0, 0.5, 0.03], facecolor=axcolor)
    rotation_speed_slider = Slider(axrotation, "Rotation Speed", -100.0, 100.0, valinit=0)
    rotation_speed_slider.on_changed(update_rotation_speed)

    
    def update(frame):
        pool.step(dt)
        balls_on_table = set(pool.balls)  # Convert the list of balls to a set for efficient lookup

        # Update positions of balls that are still on the table
        for ball, circle in ball_to_circle.items():
            if ball in balls_on_table:
                circle.center = (ball.x, ball.y)
                #ball.update_marker()  # Mise à jour de la position du marqueur
                #angle_rad = math.radians(ball.marker_angle)
                #marker_x = ball.x + math.cos(angle_rad) * ball.radius
                #marker_y = ball.y + math.sin(angle_rad) * ball.radius
                #ax.plot(marker_x, marker_y, 'o', color='black', markersize=2)

        for hole, circle in hole_circles.items():
            circle.center = (hole.x, hole.y)

        # Remove circles for balls that are no longer on the table
        for ball in list(ball_to_circle.keys()):
            if ball not in balls_on_table:
                ball_to_circle[ball].remove()  # Remove the circle from the plot
                del ball_to_circle[ball]

        if pool.white_ball not in ball_to_circle and pool.all_balls_stopped():
            circle = plt.Circle((pool.white_ball.x, pool.white_ball.y), pool.white_ball.radius, fc=pool.white_ball.color)
            ax.add_patch(circle)
            ball_to_circle[pool.white_ball] = circle

        return ball_to_circle.values()
    ani = FuncAnimation(fig, update, blit=True, interval=50, cache_frame_data=False)
    plt.show()