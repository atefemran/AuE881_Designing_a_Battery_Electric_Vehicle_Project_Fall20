# AuE881_Designing_a_Battery_Electric_Vehicle_Project
As part of AuE881 Automotive Systems Integration course, the final project was exploring the design of a battery electric vehicle using systems design approach.

## [A] Problem Statement
In this project, “design” an entire battery-electric vehicle! Given some system-level requirements, you will be asked to make design choices about six different subsystems and/or perspectives:
1. body-in-white,
2. packaging,
3. vehicle dynamics,
4. powertrain,
5. human factors,
6. and systems integration.
The overarching goal is to satisfy all the system-level requirements and ultimately, to maximize the profit for your company.

### Design problem
Design a battery-electric vehicle that maximizes profit for the enterprise and satisfies the requirements listed below. Only the details of the design as characterized by the list of design variables below need to be specified. Please include the values for these design variables for the feasible alternative(s) that you consider for your final report.

### List of design variables
* Battery capacity in [kWh]
* Choice of electric motor & inverter pair: A, B, or C (see below for motor characteristics)
* Motor and driveline configuration: FWD 1 motor, FWD 2 motors, RWD 1 motor, RWD 2 motors, AWD 2 motors, AWD 4 motors
* Final drive ratio
* 2-D vehicle layout (multiple dimensions as needed) including, L104 in [mm], L105 in [mm], L108 in [mm], H101 in [mm], H5-1 in [mm], H5-2 in [mm]
* Frame rail: height in [mm], width in [mm], and thickness in [mm]
* Suspension spring constant in [N/m]
* Suspension damping coefficient in [Ns/m]

### List of requirements:
* The vehicle shall be a medium car. Any assumptions about the packaging must take the approximate reference dimensions mentioned in the H-point book into account (see
Appendix C)
* The wheelbase, L101, shall be 2700mm (this is an artificial constraint we have imposed purely to simplify the problem).
* For the relevant packaging requirements, refer to the ranges specified as best practices in the lecture notes. This includes limits on ground lines, vision cone, head clearance, ball of foot distance, etc.
* The batteries shall be no closer than 100mm to any passenger.
* In side-view projection, passengers shall not overlap the motor(s) or transmission(s) (but ARE allowed to overlap with the driveline).
* In side-view projection, the batteries cannot overlap with the frame rail (i.e., the batteries cannot be mounted in between the frame rails).
* The vehicle shall accommodate four AM95 passengers (including the driver).
* The vehicle shall accommodate a cargo volume of 0.30 m3 and 28 kg (7kg/passenger)
* The vehicle shall be able to ascend a 12% slope at 30 km/h.
* The material stress in the frame rail shall be less than 215 MPa under normal loading conditions.
* The maximum deflection of the frame rail shall be less than 15mm.
* The vehicle shall be able to turn on a 100m diameter skidpad at up to 0.8g lateral acceleration
* To avoid the motors and inverters from burning out, the vehicle shall not exceed continuous torque limits indicated in the motor efficiency maps. We assume that continuous torque limits are at 70% of peak torque.
* Assuming the damping is at 30% of critical damping, the ride frequency (the resonant frequency) shall be equal to 1.2 Hz in the front and 1.44 Hz in the rear. Hint: These requirements depend on the weight and weight distribution of the vehicle.
* The vision cone for the driver extends no less than 6 deg down and 15 deg up relative to horizontal.
* The ground clearance of a fully rated vehicle (at GVWR) shall be no less than 130 mm.
* The angle of approach shall be no less than 12 deg.
* The angle of departure shall be no less than 12 deg.
* The ramp breakover angle shall be no less than 10 deg.
* The center of gravity should fall between 45 to 55 percent of the wheelbase measured from the front wheel axle along the vehicle length

(For the rest of the requirments and design constraints, please refer to the problem statement file)

## [B] Design Exploration
### Computational Order
The following computational order is considered:
1.	Vehicle mass
2.	Frame Rail Dimension (Max Stress, and Max deflection)
3.	CG calculations (CG to front axle over the wheelbase ratio)
4.	Vehicle costing
5.	Vehicle Performance Model (max speed, 0-100km/h time, range, 12% slope)
6.	Profit calculations
7.	Vehicle Dynamics Model (Understeering, Kus value, turning radius, and max lateral acceleration)
8.	Suspension Spring Constant and damping coefficient

### Code structure:


### Models Samples
#### The Overall Model

#### Road Load Submodel

#### Tires Reactive Forces Submodel

#### Driveline Submodel

#### Electrical Motor Submodel

#### Battery Range Calculations

More details, Please refer to the report pdf.
