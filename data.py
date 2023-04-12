from conversion import conversion

# differential ratio
differential_ratio = 4

# special number
special_num = 5252

#aerodynamic coefficient 
aerodynamic_coefficient = 0.8

# atomospheric density 
atomospheric_density = 1.269

# gravity acceleration
gravity_acceleration = 9.8

# friction coefficient
friction_coefficients = {"static": 0.6, "slide": 0.2}

# drag coefficient
drag_coefficients = {"sedan": 0.3, "sports car": 0.25, "muscle car": 0.45}

# gear ratio
gear_ratio = {"gear_1": 3, "gear_2": 2, "gear_3": 1.5, "gear_4": 1, "gear_5": 0.75, "gear_6": 0.5}

# car types
Toyota_Corolla_LE_2021 = {  "engine_power": 139,        #(hp) @ 6100rpm
                            "engine_speed": 6100, 
                            "wheel_diameter": 24.88,    #(in) P205/55HR16   
                            "mass": 2943,               #(lb)
                            "width": 70.1,              #(in)
                            "height": 56.5,             #(in)
                        }

Dodge_Challenger_1970 = {   "engine_power": 456,        #(hp) @ 6100rpm
                            "engine_speed": 5500, 
                            "wheel_diameter": 26.10,    #(in) P205/55HR16   
                            "mass": 3801,               #(lb)
                            "width": 76.5,              #(in)
                            "height": 51.0,             #(in)
                        }

Ferrari_Roma_2021     = {   "engine_power": 612,        #(hp) @ 6100rpm
                            "engine_speed": 7000, 
                            "wheel_diameter": 20,    #(in) P205/55HR16   
                            "mass": 3461,               #(lb)
                            "width": 77.7,              #(in)
                            "height": 51.2,             #(in)
                        }
