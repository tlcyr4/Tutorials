# Procedure.py
# 
# Created:  Mar 2016, M. Vegh
# Modified: 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import SUAVE
from SUAVE.Core import Units, Data
import numpy as np
import copy
from scipy import integrate
from SUAVE.Analyses.Process import Process
from SUAVE.Methods.Propulsion.turbofan_sizing import turbofan_sizing
from SUAVE.Methods.Geometry.Two_Dimensional.Cross_Section.Propulsion.compute_turbofan_geometry import compute_turbofan_geometry
#from SUAVE.Methods.Center_of_Gravity.compute_component_centers_of_gravity import compute_component_centers_of_gravity
#from SUAVE.Methods.Center_of_Gravity.compute_aircraft_center_of_gravity import compute_aircraft_center_of_gravity
from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Lift.compute_max_lift_coeff import compute_max_lift_coeff
from SUAVE.Optimization.write_optimization_outputs import write_optimization_outputs

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------   

def setup():
    
    # ------------------------------------------------------------------
    #   Analysis Procedure
    # ------------------------------------------------------------------ 
    
    # size the base config
    procedure = Process()
    procedure.simple_sizing = simple_sizing
    
    # find the weights
    procedure.weights = weight
    # finalizes the data dependencies
    procedure.finalize = finalize
    
    # performance studies
    procedure.missions                   = Process()
    procedure.missions.design_mission    = design_mission

    # post process the results
    procedure.post_process = post_process
        
    # done!
    return procedure


# ----------------------------------------------------------------------        
#   Target Range Function
# ----------------------------------------------------------------------    

def find_target_range(nexus,mission):
    
    segments=mission.segments
    cruise_altitude=mission.segments['climb_5'].altitude_end
    climb_1=segments['climb_1']
    climb_2=segments['climb_2']
    climb_3=segments['climb_3']
    climb_4=segments['climb_4']
    climb_5=segments['climb_5']
  
    descent_1=segments['descent_1']
    descent_2=segments['descent_2']
    descent_3=segments['descent_3']

    x_climb_1=climb_1.altitude_end/np.tan(np.arcsin(climb_1.climb_rate/climb_1.air_speed))
    x_climb_2=(climb_2.altitude_end-climb_1.altitude_end)/np.tan(np.arcsin(climb_2.climb_rate/climb_2.air_speed))
    x_climb_3=(climb_3.altitude_end-climb_2.altitude_end)/np.tan(np.arcsin(climb_3.climb_rate/climb_3.air_speed))
    x_climb_4=(climb_4.altitude_end-climb_3.altitude_end)/np.tan(np.arcsin(climb_4.climb_rate/climb_4.air_speed))
    x_climb_5=(climb_5.altitude_end-climb_4.altitude_end)/np.tan(np.arcsin(climb_5.climb_rate/climb_5.air_speed))
    x_descent_1=(climb_5.altitude_end-descent_1.altitude_end)/np.tan(np.arcsin(descent_1.descent_rate/descent_1.air_speed))
    x_descent_2=(descent_1.altitude_end-descent_2.altitude_end)/np.tan(np.arcsin(descent_2.descent_rate/descent_2.air_speed))
    x_descent_3=(descent_2.altitude_end-descent_3.altitude_end)/np.tan(np.arcsin(descent_3.descent_rate/descent_3.air_speed))
    
    cruise_range=mission.design_range-(x_climb_1+x_climb_2+x_climb_3+x_climb_4+x_climb_5+x_descent_1+x_descent_2+x_descent_3)
  
    segments['cruise'].distance=cruise_range
    
    return nexus

# ----------------------------------------------------------------------        
#   Design Mission
# ----------------------------------------------------------------------    
def design_mission(nexus):
    
    mission = nexus.missions.base
    mission.design_range = 1500.*Units.nautical_miles
    find_target_range(nexus,mission)
    results = nexus.results
    results.base = mission.evaluate()
    
    return nexus



# ----------------------------------------------------------------------        
#   Sizing
# ----------------------------------------------------------------------    

def simple_sizing(nexus):
    configs=nexus.vehicle_configurations
    base=configs.base
    
    for config in configs:   
        if config.tag == 'initial': # ?????
            continue  
        for wing in config.wings:
            
            # This procedure is created under the assumption that only the main wing is being modified
            if wing.tag != 'main_wing':
                #wing = SUAVE.Methods.Geometry.Two_Dimensional.Planform.wing_planform(wing)
                if wing.tag == 'horizontal_stabilizer':
                    wing.areas.reference        = configs.initial.wings[wing.tag].areas.reference*config.wings['main_wing'].areas.reference/configs.initial.wings['main_wing'].areas.reference                
                wing.areas.wetted   = 2.0 * wing.areas.reference
                wing.areas.exposed  = 0.8 * wing.areas.wetted
                wing.areas.affected = 0.6 * wing.areas.reference                
                continue
            
            ## For testing OpenVSP issue
            wing.areas.wetted   = 2.0 * wing.areas.reference

            ## For testing OpenVSP issue
            wing.areas.wetted   = 2.0 * wing.areas.reference
            #wing.areas.exposed  = 0.8 * wing.areas.wetted
            #wing.areas.affected = 0.6 * wing.areas.reference                
            #continue            

            # *** Initial values of wing geometry ***
            
            initial_origin = configs.initial.wings[wing.tag].origin
       
            # initial reference area 
            S0   = configs.initial.wings[wing.tag].areas.reference
            # full span 
            b0_full_span  = configs.initial.wings[wing.tag].spans.projected 
             
            # intial chord,  relative chond sprd orgin aan of each segment of the wing stored in array 
            r0 = []  #r[0] is root chord, r[i=1,2..n] = end chord of segment i
            b0 = []  #b[0] is semispan, b[i=1,2..n] = span of segment i
            rs0 = [] #rs0[0] = 0, rs0[i] = x distance from segment root origin to segment tip origin
            n_segments = len(wing.Segments.keys())
            #print n_segments
            for i_segs in xrange(n_segments):
                if i_segs == 0:
                    r0.append(configs.initial.wings[wing.tag].chords.root)
                    b0.append(configs.initial.wings[wing.tag].spans.projected*0.5)
                    rs0.append(0)
                else:
                    b0.append((configs.initial.wings[wing.tag].Segments[i_segs].percent_span_location - configs.initial.wings[wing.tag].Segments[i_segs-1].percent_span_location)*b0[0])
                    r0.append(r0[0]*configs.initial.wings[wing.tag].Segments[i_segs].root_chord_percent)
                    rs0.append(b0[i_segs]*np.tan(configs.initial.wings[wing.tag].Segments[i_segs-1].sweeps.quarter_chord))
              
            anchor_x = initial_origin[0] + .25*r0[0]
       

            # *** modified geometry of wing ***            
            # unpack
            if wing.tag == 'horizontal_stabilizer':
                sref        = configs.initial.wings[wing.tag].areas.reference*config.wings['main_wing'].areas.reference/configs.initial.wings['main_wing'].areas.reference
            sref        = wing.areas.reference
            taper       = wing.taper
            sweep       = wing.sweeps.quarter_chord
            ar          = wing.aspect_ratio
            t_c_w       = wing.thickness_to_chord
            dihedral    = wing.dihedral 
            vertical    = wing.vertical
            symmetric   = wing.symmetric
            origin      = wing.origin
            
            # check 
            wing_tip_constraint = 0.

            # calculate
            span       = (ar*sref)**.5
            chord_root = r0[0]*sref/S0*b0_full_span /span
            chord_tip  = 0.0 #imposed condition for chord tip 
            b_full_span = span 
            
            r = [] #r[0] is root chord, r[i] = end chord of segment i
            b = [] #b[0] is semispan, b[i] = span of segment i
            for i_segs in  xrange(n_segments):
                if i_segs == 0: 
                    r.append(chord_root)
                    b.append(0.5*span)
                elif i_segs == n_segments: 
                    r.append(chord_tip)
                else:
                    r.append(r0[i_segs]*(r[0]/r0[0]))
                    b.append(b0[i_segs]*(b_full_span/b0_full_span ))
                
            trap_area = 0
            for i_segs in xrange(n_segments-1):
                trap_area = trap_area + ((r[i_segs] + r[i_segs+1])/2)*b[i_segs+1]


            # new wetted area 
            swet = 2.*2.*trap_area*(1.0 + 0.2*t_c_w) # top/bottom + left/right 
            
            # new mean aerodynamic chord  ????
            mac_arg = 0
            for i_segs in xrange(n_segments-1):
                func = lambda y: (r[i_segs] - ((r[i_segs] - r[i_segs+1])/b[i_segs+1])*y)**2
                mac_arg =  mac_arg + integrate.quad(func, 0, b[i_segs+1])[0] 
            
            mac = mac_arg*(2/sref)
            
            rs = []
            theta = []
            # calculate sweeps
            for i_segs in xrange(n_segments):
                rs.append(rs0[i_segs]*(r[0]/r0[0]))
            
            for i_segs in xrange(n_segments-1):           
                theta.append(np.arctan(rs[i_segs+1]/b[i_segs+1]))
            
            theta.append(0)
            
            # estimating aerodynamic center coordinates
            yb = []
            A = []
            for i_segs in xrange(n_segments-1):           
                yb.append( (r[i_segs]+2*r[i_segs+1])/(3*(r[i_segs] + r[i_segs + 1]))*b[i_segs+1])
                A.append((r[i_segs] + r[i_segs+1])/2.*b[i_segs+1])

            y_coord = (np.dot(yb,A))/sref*2.  #??????
            x_coord = y_coord*np.tan(theta[0]) #?????
            z_coord = y_coord * np.tan(dihedral) # does not consider possibility of multiple sections
                
            if vertical:
                temp    = y_coord * 1.
                y_coord = z_coord * 1.
                z_coord = temp
        
            if symmetric:
                y_coord = 0.   
                
          
            # update
            wing.chords.root                = chord_root
            wing.total_length               = chord_root
            wing.chords.tip                 = chord_tip
            wing.chords.mean_aerodynamic    = mac
            wing.areas.wetted               = swet
            wing.spans.projected            = span
            wing.aerodynamic_center         = [x_coord , y_coord, z_coord]    
            
            for i_segments in xrange(n_segments):
                wing.Segments[i_segs].sweeps.quarter_chord = theta[i_segs]
                wing.Segments[i_segs].root_chord_percent    = r[i_segs]/r[0]   
                        
            wing.origin[0] = anchor_x - .25*r[0]
                
            wing.areas.exposed  = 0.8 * wing.areas.wetted
            wing.areas.affected = 0.6 * wing.areas.reference  
    
    #find conditions
    air_speed   = nexus.missions.base.segments['cruise'].air_speed 
    altitude    = nexus.missions.base.segments['climb_5'].altitude_end
    atmosphere  = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    
    freestream  = atmosphere.compute_values(altitude)
    freestream0 = atmosphere.compute_values(6000.*Units.ft)  #cabin altitude
    
    
    diff_pressure         = np.max(freestream0.pressure-freestream.pressure,0)
    fuselage              = base.fuselages['fuselage']
    fuselage.differential_pressure = diff_pressure 
    
    #now size engine
    mach_number        = air_speed/freestream.speed_of_sound
    
    #now add to freestream data object
    freestream.velocity    = air_speed
    freestream.mach_number = mach_number
    freestream.gravity     = 9.81
    
    conditions             = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()   #assign conditions in form for propulsor sizing
    conditions.freestream  = freestream
    
    
    for config in configs:
        config.wings.horizontal_stabilizer.areas.reference = (26.0/92.0)*config.wings.main_wing.areas.reference
            
        #for wing in config.wings:
            
            #wing = SUAVE.Methods.Geometry.Two_Dimensional.Planform.wing_planform(wing)
            
            #wing.areas.exposed  = 0.8 * wing.areas.wetted
            #wing.areas.affected = 0.6 * wing.areas.reference
            


        fuselage              = config.fuselages['fuselage']
        fuselage.differential_pressure = diff_pressure 
        
        turbofan_sizing(config.propulsors['turbofan'], mach_number = mach_number, altitude = altitude)
        compute_turbofan_geometry(config.propulsors['turbofan'], conditions)


    ## ------------------------------------------------------------------
    ##   Landing Configuration
    ## ------------------------------------------------------------------
    #landing = nexus.vehicle_configurations.landing
    #landing_conditions = Data()
    #landing_conditions.freestream = Data()

    ## landing weight
    #landing.mass_properties.landing = 0.85 * config.mass_properties.takeoff
    
    ## Landing CL_max
    #altitude = nexus.missions.base.segments[-1].altitude_end
    #atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    #freestream_landing  = atmosphere.compute_values(0.)
    ##p, T, rho, a, mu = atmosphere.compute_values(0.)
    #landing_conditions.freestream.velocity           = nexus.missions.base.segments['descent_3'].air_speed
    #landing_conditions.freestream.density            = freestream_landing.density
    #landing_conditions.freestream.dynamic_viscosity  = freestream_landing.dynamic_viscosity
    #CL_max_landing,CDi = compute_max_lift_coeff(landing,landing_conditions)
    #landing.maximum_lift_coefficient = CL_max_landing
    
    ##Takeoff CL_max
    #takeoff = nexus.vehicle_configurations.takeoff
    #takeoff_conditions = Data()
    #takeoff_conditions.freestream = Data()    
    #altitude = nexus.missions.base.airport.altitude
    #freestream_takeoff  = atmosphere.compute_values(altitude)
   
    ##p, T, rho, a, mu = atmosphere.compute_values(altitude)
    #takeoff_conditions.freestream.velocity           = nexus.missions.base.segments.climb_1.air_speed
    #takeoff_conditions.freestream.density            = freestream_takeoff.density
    #takeoff_conditions.freestream.dynamic_viscosity  = freestream_takeoff.dynamic_viscosity 
    #max_CL_takeoff,CDi = compute_max_lift_coeff(takeoff,takeoff_conditions) 
    #takeoff.maximum_lift_coefficient = max_CL_takeoff
    
    ##Base config CL_max
    #base = nexus.vehicle_configurations.base
    #base_conditions = Data()
    #base_conditions.freestream = takeoff_conditions.freestream   
    #max_CL_base,CDi = compute_max_lift_coeff(base,base_conditions) 
    #base.maximum_lift_coefficient = max_CL_base    
    ## done!
    
    return nexus

# ----------------------------------------------------------------------        
#   Weights
# ----------------------------------------------------------------------    

def weight(nexus):
    vehicle=nexus.vehicle_configurations.base

    # weight analysis
    weights = nexus.analyses.base.weights.evaluate()
   
    '''
    compute_component_centers_of_gravity(vehicle)
    nose_load_fraction=.06
    compute_aircraft_center_of_gravity(vehicle,nose_load_fraction)
    '''
    
    weights = nexus.analyses.cruise.weights.evaluate()
    vehicle.mass_properties.breakdown = weights
    weights = nexus.analyses.landing.weights.evaluate()
    weights = nexus.analyses.takeoff.weights.evaluate()
    weights = nexus.analyses.short_field_takeoff.weights.evaluate()
    
    empty_weight    =vehicle.mass_properties.operating_empty
    passenger_weight=vehicle.passenger_weights.mass_properties.mass 
    for config in nexus.vehicle_configurations:
        #config.mass_properties.max_zero_fuel                = empty_weight+passenger_weight
        config.mass_properties.zero_fuel_center_of_gravity  = vehicle.mass_properties.zero_fuel_center_of_gravity
        config.fuel                                         = vehicle.fuel
       
    return nexus


# ----------------------------------------------------------------------
#   Finalizing Function (make part of optimization nexus)[needs to come after simple sizing doh]
# ----------------------------------------------------------------------    

def finalize(nexus):
    
    nexus.analyses.finalize()   
    
    return nexus         


    
# ----------------------------------------------------------------------
#   Post Process Results to give back to the optimizer
# ----------------------------------------------------------------------   

def post_process(nexus):
    
    # Unpack data
    vehicle                           = nexus.vehicle_configurations.base
    
    '''
    print 'base.mass_properties.takeoff = ', vehicle.mass_properties.takeoff
    print 'takeoff.mass_properties.takeoff = ',  nexus.vehicle_configurations.takeoff.mass_properties.takeoff
    print 'vehicle.mass_properties.empty = ', vehicle.mass_properties.operating_empty
    '''
    results                           = nexus.results
    summary                           = nexus.summary
    missions                          = nexus.missions  
    nexus.total_number_of_iterations +=1
    # Static stability calculations
    CMA = -10.
    for segment in results.base.segments.values():
        max_CMA=np.max(segment.conditions.stability.static.cm_alpha[:,0])
        if max_CMA>CMA:
            CMA=max_CMA
            

            
    summary.static_stability = CMA
    
    #throttle in design mission
    max_throttle=0
    for segment in results.base.segments.values():
        max_segment_throttle = np.max(segment.conditions.propulsion.throttle[:,0])
        if max_segment_throttle > max_throttle:
            max_throttle = max_segment_throttle

            
    summary.max_throttle = max_throttle
    
    # Fuel margin and base fuel calculations
    operating_empty          = vehicle.mass_properties.operating_empty
    payload                  = vehicle.passenger_weights.mass_properties.mass 
    design_landing_weight    = results.base.segments[-1].conditions.weights.total_mass[-1]
    design_takeoff_weight    = vehicle.mass_properties.takeoff
    max_takeoff_weight       = nexus.vehicle_configurations.takeoff.mass_properties.max_takeoff
    zero_fuel_weight         = payload+operating_empty
    
    summary.max_zero_fuel_margin    = (design_landing_weight - zero_fuel_weight)/zero_fuel_weight
    summary.base_mission_fuelburn   = design_takeoff_weight - results.base.segments['descent_3'].conditions.weights.total_mass[-1]
 
  

    hf = vehicle.fuselages.fuselage.heights.at_wing_root_quarter_chord
    wf = vehicle.fuselages.fuselage.width
    Lf = vehicle.fuselages.fuselage.lengths.total
    Sw = vehicle.wings.main_wing.areas.reference
    cw = vehicle.wings.main_wing.chords.mean_aerodynamic
    b  = vehicle.wings.main_wing.spans.projected
    Sh = vehicle.wings.horizontal_stabilizer.areas.reference
    Sv = vehicle.wings.vertical_stabilizer.areas.reference
    lh = vehicle.wings.horizontal_stabilizer.origin[0] + vehicle.wings.horizontal_stabilizer.aerodynamic_center[0] - vehicle.mass_properties.center_of_gravity[0]
    lv = vehicle.wings.vertical_stabilizer.origin[0] + vehicle.wings.vertical_stabilizer.aerodynamic_center[0] - vehicle.mass_properties.center_of_gravity[0]

    
    #when you run want to output results to a file
    filename = 'results.txt'
    write_optimization_outputs(nexus, filename)
    '''
    unscaled_inputs = nexus.optimization_problem.inputs[:,1] #use optimization problem inputs here
    input_scaling   = nexus.optimization_problem.inputs[:,3]
    scaled_inputs   = unscaled_inputs/input_scaling
    problem_inputs=[]
    
    for value in unscaled_inputs:
        problem_inputs.append(value) 
    file=open('results.txt' , 'ab')
    file.write('iteration = ')
    file.write(str(nexus.iteration_number))
    file.write('fuel weight = ')
    file.write(str( summary.base_mission_fuelburn))
  
    file.write(', inputs = ')
    file.write(str(problem_inputs))
    
    file.write('\n') 
    file.close()
    '''
    
    
    
    
    return nexus    
