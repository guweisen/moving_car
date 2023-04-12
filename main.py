import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

import math

import data
from conversion import conversion

# Checking the velocities and net force to see if it makes sense


# ─── CONSTANTS ──────────────────────────────────────────────────────────────────
gear_ratio = data.gear_ratio
differential_ratio = data.differential_ratio
fifty_two_fifty_two = data.special_num
aerodynamic_coefficient = data.aerodynamic_coefficient
atomospheric_density = data.atomospheric_density
drag_coefficients = data.drag_coefficients
gravity_acceleration = data.gravity_acceleration
friction_coefficients = data.friction_coefficients

#
# ──────────────────────────────────────────────────────────────────────────────── I ──────────
#   :::::: C A R   I N F O   C A L C U L A T I O N : :  :   :    :     :        :          :
# ──────────────────────────────────────────────────────────────────────────────────────────
#

# ─── CAR INFO ───────────────────────────────────────────────────────────────────
class Car:
    def __init__(self, user_input):
        self.user_input = user_input  
        self.fixed_values = {}

        self.gear_top_speed_mps_dict = {}       #(m/s)
        self.gear_acceleration_mps2_dict = {}   #(m/s^2)

        self.gear_wheel_torque_dict = {}        #(Nm)
        self.gear_driving_force_dict = {}       #(N)

        self.mass = user_input["mass"]                      #(lb)
        self.engine_power = user_input["engine_power"]      #(hp)
        self.engine_speed = user_input["engine_speed"]      #(rpm)
        self.wheel_diameter = user_input["wheel_diameter"]  #(in)
 
# ─── FIXED VALUES CALCULATION ──────────────────────────────────────────────────────
    def fixed_values_calculation(self):
        width_height = user_input["width"] * user_input["height"]
        frontal_area = (width_height * aerodynamic_coefficient) * conversion["in to m"]    #(m) 
        self.fixed_values["frontal_area"] = round(frontal_area, 4)

        wheel_circumference = self.wheel_diameter * math.pi                 #(inch) distance travel every rev
        self.fixed_values["wheel_circumference"] = round(wheel_circumference, 4)

        engine_speed_hr = self.engine_speed * 60                            #(rph)rev per hour 
        self.fixed_values["engine_speed_hr"] = engine_speed_hr

        engine_torque = (self.engine_power * fifty_two_fifty_two / self.engine_speed)   
        self.fixed_values["engine_torque"] = round(engine_torque, 4)        #(lb-ft)

        normal_force = (self.mass * conversion["lb to kg"]) * gravity_acceleration      
        self.fixed_values["normal_force"] = round(normal_force, 4)          #(N)
 
        friction_force_static =  friction_coefficients["static"] * normal_force 
        friction_force_slide =  friction_coefficients["slide"] * normal_force 
        self.fixed_values["friction_force_static"] = round(friction_force_static, 4)    #(N)
        self.fixed_values["friction_force_slide"] = round(friction_force_slide, 4)      #(N)

# ─── GEAR CALCULATIONS ──────────────────────────────────────────────────────────
    def gear_calculation(self):
        for gears in gear_ratio:
            gear = gear_ratio[gears]                                # gear ratio of a specific gear
            gear_differential_ratio = differential_ratio * gear     # calculate the gear and differential ratio

            gear_top_speed_inch_hr = (self.fixed_values["engine_speed_hr"] / gear_differential_ratio)* self.fixed_values["wheel_circumference"]
            gear_top_speed_mps = round(gear_top_speed_inch_hr * conversion["in/hr to m/s"], 4)          #(m/s)
            self.gear_top_speed_mps_dict[gears+'_top_speed'] = gear_top_speed_mps

            wheel_torque = (self.fixed_values["engine_torque"] * gear_differential_ratio) / ((self.wheel_diameter / 2) / 12)
            wheel_torque = round(wheel_torque * conversion["lb-ft to Nm"], 2)       #(lb-ft)
            self.gear_wheel_torque_dict[gears+'_wheel_torque'] = wheel_torque               

            wheel_force = wheel_torque / ((self.wheel_diameter / 2) / 12) 
            wheel_force = round(wheel_force, 4)                                     #(lbf)
            driving_force = round(wheel_force * conversion["lbf to N"], 4)          #(N)
            self.gear_driving_force_dict[gears+"_driving_force"] = driving_force
        
            # acceleration = wheel_force / self.mass                                  #(g's)
            # acceleration_mps2 = round(acceleration * conversion["g's to m/s^2"], 4)
            # self.gear_acceleration_mps2_dict[gears+"_acceleration"] = acceleration_mps2 
    
    def car_calculation(self):
        self.fixed_values_calculation()
        self.gear_calculation()

# ────────────────────────────────────────────────────────────────────────────────

# ─── SIMULATION ─────────────────────────────────────────────────────────────────
class Simulation:
    def __init__(self, car, gas_pedal_pressed):
        self.car = car
        self.gas_pedal_pressed = gas_pedal_pressed

        self.current_velocity = 0 
        self.previous_velocity = 0
        self.acceleration = 0
        self.terminal_velocity = False
        
        self.air_resistance = 0
        self.total_friction = 0
        self.net_force = 0
        self.change_in_momentum_force = 0

        self.start_gear = 1
        self.gear_number = "gear_" + str(self.start_gear) + "_driving_force"

        self.current_time = 0
        self.time_step = 0.1

    def gear_change_check(self, car):
        gear_top_speed_number = "gear_" + str(self.start_gear) + "_top_speed"
        if self.previous_velocity >= self.car.gear_top_speed_mps_dict[gear_top_speed_number] and self.start_gear < 6:
            self.start_gear += 1
        elif self.previous_velocity < self.car.gear_top_speed_mps_dict[gear_top_speed_number] and self.start_gear > 1:
            self.start_gear -= 1

    def net_force_calculation(self, car):
        self.air_resistance = (1/2) * atomospheric_density * drag_coefficient * (self.previous_velocity ** 2) * self.car.fixed_values["frontal_area"]
        self.total_friction = self.air_resistance + self.car.fixed_values["friction_force_slide"]

        if self.gas_pedal_pressed:
            if self.current_velocity == 0:
                self.net_force = self.car.gear_driving_force_dict[self.gear_number] - self.car.fixed_values["friction_force_static"]    # gas pedal pressed and car hasn't moved yet
            else:
                self.net_force = (self.car.gear_driving_force_dict[self.gear_number] + self.change_in_momentum_force) - self.total_friction     # gas pedal pressed and car has been moving
        else:
            self.net_force = self.change_in_momentum_force - self.total_friction    # gas pedal not pressed so no force from the engine
    
    def terminal_velocity_check(self):
        if self.gas_pedal_pressed and self.net_force <= self.total_friction and self.acceleration == 0:
            self.terminal_velocity = True
        else:
            self.terminal_velocity = False
            
    def update_gas(self, car):
        self.gear_change_check(car)
        self.net_force_calculation(car)
        self.terminal_velocity_check()
        
        if self.current_velocity == 0:             # initial state
            self.acceleration = self.net_force / self.car.mass
        elif self.net_force > self.total_friction:    # gas pedal pressed and car is moving with +a
            self.acceleration = self.net_force / self.car.mass
        elif self.gas_pedal_pressed == True:
            self.acceleration = 0                   # gas pedal pressed but reached terminal velocity

        self.previous_velocity = self.current_velocity
        if self.previous_velocity < self.car.gear_top_speed_mps_dict["gear_6_top_speed"]:
            self.current_velocity += self.acceleration * self.time_step
        else:
            self.previous_velocity = self.current_velocity

        self.change_in_momentum_force = (self.car.mass * (self.current_velocity - self.previous_velocity)) / self.time_step
    
        self.current_time += self.time_step 

    def update_no_gas(self, car):
        self.gear_change_check(car)
        self.net_force_calculation(car)
        
        self.acceleration = self.net_force / self.car.mass   # gas pedal is not pressed, acceleration will decrease
        
        self.previous_velocity = self.current_velocity
        self.current_velocity += self.acceleration * self.time_step

        self.change_in_momentum_force = (self.car.mass * (self.current_velocity - self.previous_velocity)) / self.time_step
    
        self.current_time += self.time_step         

    def print_info(self, car) :
        print(f"""
        UPDATE NO GAS
        Gear: {self.start_gear}
        Gas Pedal: {self.gas_pedal_pressed}
        Current Time: {self.current_time}
        Air Resistance: {self.air_resistance}
        Change in Momentum Force: {self.change_in_momentum_force}
        Net Force: {self.net_force}
        Acceleration: {self.acceleration}
        Previous Velocity: {self.previous_velocity}
        Current Velocity: {self.current_velocity}
        Terminal Velocity: {self.terminal_velocity}
        """)

user_input = input("----- Enter a letter for the car type -----\n a) Toyota Corolla 2021 LE \n b) Dodge Challenger 1970 \n c) Ferrari Roma 2021")
user_input = user_input.lower()

if user_input == "a":
    user_input = data.Toyota_Corolla_LE_2021
    drag_coefficient = drag_coefficients["sedan"]
elif user_input == "b":
    user_input = data.Dodge_Challenger_1970
    drag_coefficient = drag_coefficients["muscle car"]
elif user_input == "c":
    user_input = data.Ferrari_Roma_2021
    drag_coefficient = drag_coefficients["sports car"]

# ─── MATPLOT VALUES ─────────────────────────────────────────────────────────────────
velocities_mps = [0]
velocities_kmph = [0]
accelerations = [0]

net_forces = [0]
total_frictions = [0]

time = [0]
terminal_velocities = []

num_steps = 50

car = Car(user_input)
car.car_calculation()

sim = Simulation(car, gas_pedal_pressed=True)


def update_data():
    net_force = sim.net_force
    net_forces.append(net_force)

    total_friction = sim.total_friction
    total_frictions.append(total_friction)

    velocity = sim.current_velocity
    velocities_mps.append(velocity)
    velocities_kmph.append(velocity * conversion["m/s to km/h"])

    acceleration = sim.acceleration
    accelerations.append(acceleration)

    velocity_time = sim.current_time
    time.append(velocity_time)
    
fig, axs = plt.subplots(3,2)

# ─── AXES CREATION ─────────────────────────────────────────────────────────────────
velocities_mps_ax = axs[0, 0]
velocities_kmph_ax = axs[1, 0]
accelerations_ax = axs[2, 0]

bigger_velocities_kmph_ax = axs[0, 1]
total_frictions_ax = axs[1, 1]
net_forces_ax = axs[2, 1]

# ─── AXES TITLES ─────────────────────────────────────────────────────────────────
velocities_mps_ax.set_title("Velocity(m/s)")
velocities_mps_ax.set_ylabel("(m/s)")
velocities_kmph_ax.set_title("Velocity(km/h)")
velocities_kmph_ax.set_ylabel("(km/h)")
accelerations_ax.set_title("Acceleration(m/s^2)")
accelerations_ax.set_ylabel("(m/s^2)")

bigger_velocities_kmph_ax.set_title("Velocity(km/h)")
bigger_velocities_kmph_ax.set_ylabel("(km/h)")

total_frictions_ax.set_title("Drag Force(N)")
total_frictions_ax.set_ylabel("(N)")

net_forces_ax.set_title("Net Force(N)")
net_forces_ax.set_ylabel("(N)")

# ─── VALUES TO TIME PLOTS ─────────────────────────────────────────────────────────────────
velocities_mps_time, = velocities_mps_ax.plot(time, velocities_mps)
velocities_kmph_time, = velocities_kmph_ax.plot(time, velocities_kmph)
accelerations_time, = accelerations_ax.plot(time, accelerations)

bigger_velocities_kmph_time, = bigger_velocities_kmph_ax.plot(time, velocities_kmph)
bigger_velocities_kmph_ax.set_xlim(0, sim.time_step* num_steps + 2)
bigger_velocities_kmph_ax.set_ylim(0, (car.gear_top_speed_mps_dict["gear_6_top_speed"] * conversion["m/s to km/h"]) * 1.5)

net_forces_time, = net_forces_ax.plot(time, net_forces)
total_frictions_time, = total_frictions_ax.plot(time, total_frictions)

# ─── ARTISTS ─────────────────────────────────────────────────────────────────
def artists(fig, axs):
    return velocities_mps_time, velocities_kmph_time, accelerations_time, bigger_velocities_kmph_time, net_forces_time, total_frictions_time

# ─── ANIMATE FUNCTION ─────────────────────────────────────────────────────────────────
def animate(frame):
    if velocities_mps[-1] < 0:
        animation.event_source.stop()
        # plt.close(fig)
    elif len(terminal_velocities) < 10:
        sim.update_gas(car)
        # sim.print_info(car)
        update_data()
        if sim.terminal_velocity:
            terminal_velocities.append(1)
    else:
        sim.gas_pedal_pressed = False
        sim.update_no_gas(car)
        update_data()

# ─── UPDATING VALUES ─────────────────────────────────────────────────────────────────
    velocities_mps_time.set_data(time, velocities_mps)
    velocities_kmph_time.set_data(time, velocities_kmph)
    accelerations_time.set_data(time, accelerations)
    
    bigger_velocities_kmph_time.set_data(time, velocities_kmph)
    net_forces_time.set_data(time, net_forces)
    total_frictions_time.set_data(time, total_frictions)

# ─── PLOTS LIMITS ─────────────────────────────────────────────────────────────────
    for ax in range(len(axs)):
        axs[ax, 0].relim()
        axs[ax, 0].autoscale_view()

    net_forces_ax.relim()
    net_forces_ax.autoscale_view()

    total_frictions_ax.relim()
    total_frictions_ax.autoscale_view()

    fig.subplots_adjust(hspace=0.5)
    fig.canvas.draw_idle()

animation = FuncAnimation(fig, animate, frames=100, interval=200, init_func=lambda: artists(fig, axs)) 
plt.show()

