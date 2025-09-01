use std::{error::Error, f64::{NAN, consts::PI}};
use std::time::Instant;
use csv::Writer;

//do not touch constants
const DIMENTIONALITY: usize = 3; //num of spatial dimensions, higher or lower dimensionality not implemented, for 2D use set z value = 0.0
type FRes = f64; //can be changed to f32, f128 etc depending on need
type FixedVec = [FRes; DIMENTIONALITY]; //vector fixed in size (hence the name)
const GRAVITATIONAL_CONST: FRes = 6.674_3e-11;
const ORBITS_TO_COMPLETE: i32 = 5; //I would allow for this to be a higher or lower number however at this point I think it's not worth the time so don't change this

//object constants
const EARTH_MASS: FRes = 5.972e24; //mass inspired by earth
const HUBBLE_MASS: FRes = 11_110_000.0; //mass inspired by hubble telescope
const RADIUS_EARTH: FRes = 6_378_000.0;
const ORBITAL_RADIUS: FRes = 100_000.0 + RADIUS_EARTH;

//relatively "adjustable" constants (adjustable before running program of course)
const TIME_STEP: FRes = 0.01;
const V_ESCAPE: FRes = 11093.2113505; //these are declared here as they may be useful, they were calculated using the formulae for their respective velocity
const V_CIRC_ORBIT: FRes = 7844.08497109;
const NUM_VELOCITIES: usize = 400;
const MAX_TAN_VELOCITY: FRes = V_CIRC_ORBIT*{5.0/4.0};
const MIN_TAN_VELOCITY: FRes = -MAX_TAN_VELOCITY;
const V_STEP: FRes = (MAX_TAN_VELOCITY-MIN_TAN_VELOCITY)/{NUM_VELOCITIES as f64};

struct GObject {
    position: FixedVec,
    velocity: FixedVec,
    acceleration: FixedVec,
    mass: FRes,
}

#[derive(serde::Serialize)]
struct CSVDataInstance {
    vel_tan: FRes,
    time_1: FRes,
    time_2: FRes,
    time_3: FRes,
    time_4: FRes,
    time_5: FRes,
}
impl CSVDataInstance{
    fn dump_time_array(&mut self, array: &[FRes; ORBITS_TO_COMPLETE as usize]){
        self.time_1 = array[0];
        self.time_2 = array[1];
        self.time_3 = array[2];
        self.time_4 = array[3];
        self.time_5 = array[4];
    }
}

fn subtract_vecs(vec1: &FixedVec, vec2: &FixedVec, vec_out: &mut FixedVec){
    vec_out[0] = vec1[0] - vec2[0];
    vec_out[1] = vec1[1] - vec2[1];
    vec_out[2] = vec1[2] - vec2[2];
}

fn add_vec_to_obj(vec1: &FixedVec, object_vec: &mut FixedVec){
    object_vec[0] += vec1[0];
    object_vec[1] += vec1[1];
    object_vec[2] += vec1[2];
}

fn set_squared_dot_product(vec1: &FixedVec, vec2: &FixedVec) -> FRes{
    let squared_dprod: FRes = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    return squared_dprod;
}

fn multiply_passed_vec(vector: &FixedVec, mult: FRes, vec_out: &mut FixedVec){
    vec_out[0] = vector[0] * mult;
    vec_out[1] = vector[1] * mult;
    vec_out[2] = vector[2] * mult;
}

fn scale_vec(vector: &mut FixedVec, mult: FRes){
    vector[0] *= mult;
    vector[1] *= mult;
    vector[2] *= mult;
}

fn transfer_to_unit(dot_product: FRes, transfer: &mut FixedVec){
    transfer[0] /= dot_product;
    transfer[1] /= dot_product;
    transfer[2] /= dot_product;
}

fn set_accelerations(position1: &FixedVec, position2: &FixedVec, 
    accleration1: &mut FixedVec, acceleration2: &mut FixedVec,
    mass1: &FRes, mass2: &FRes, 
    transfer_vector /*used for transferring vectors to objects*/: &mut FixedVec){
    subtract_vecs( position1, position2, transfer_vector); //transfer acts as displacement vector
    let squared_dot: FRes = set_squared_dot_product(transfer_vector, transfer_vector); //pass in reference instead of allocating each time?
    transfer_to_unit(squared_dot.sqrt(), transfer_vector); //transfer vector set to unit vector (might cause issues?)
    multiply_passed_vec(transfer_vector, -GRAVITATIONAL_CONST*(mass2)/squared_dot, accleration1);
    scale_vec(transfer_vector, -1.0); //flip for other object
    multiply_passed_vec(transfer_vector, -GRAVITATIONAL_CONST*(mass1)/squared_dot, acceleration2);

}

fn leap_frog(object1: &mut GObject, object2: &mut GObject, transfer_vector /*used for transferring vectors to objects*/: &mut FixedVec){
    //setting v_half for objects eg. v_{i+1/2} = v_i + 1/2 * a_i * delta_t
    multiply_passed_vec(&object1.acceleration, 0.5*TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object1.velocity);
    multiply_passed_vec(&object2.acceleration, 0.5*TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object2.velocity);

    //updating positions: x_{i+1} = x_i + v_{i+1/2} * delta_t
    multiply_passed_vec(&object1.velocity, TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object1.position);
    multiply_passed_vec(&object2.velocity, TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object2.position);

    //recalculate accelerations
    set_accelerations(&object1.position, &object2.position, 
        &mut object1.acceleration, &mut object2.acceleration, 
        &object1.mass, &object2.mass, transfer_vector);
    
    //set v_{i+1} = v_{i+1/2} + 1/2 * a_{i+1} * delta_t
    multiply_passed_vec(&object1.acceleration, 0.5*TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object1.velocity);
    multiply_passed_vec(&object2.acceleration, 0.5*TIME_STEP, transfer_vector);
    add_vec_to_obj(transfer_vector, &mut object2.velocity);
}

fn get_delta_theta(old_displacement: &mut FixedVec, new_displacement: &mut FixedVec) -> FRes{
    let mut delta_theta: FRes = set_squared_dot_product(old_displacement, new_displacement);
    delta_theta = delta_theta/(set_squared_dot_product(old_displacement, old_displacement).sqrt()*set_squared_dot_product(new_displacement, new_displacement).sqrt());
    return delta_theta.acos();
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut csv_writer = Writer::from_path("orbital_periods_output.csv")?;
    let mut earth = GObject { 
        position: [0.0, 0.0, 0.0], 
        velocity: [0.0, 0.0, 0.0],
        acceleration: [0.0, 0.0, 0.0],
        mass: EARTH_MASS
    };
    let earth_ref: &mut GObject = &mut earth;

    let mut satellite_object= GObject{
        position: [0.0, 100_000.0+RADIUS_EARTH, 0.0],
        velocity: [0.0, 0.0, 0.0],
        acceleration: [0.0, 0.0, 0.0],
        mass: HUBBLE_MASS
    };
    let satellite_ref: &mut GObject = &mut satellite_object;
    let transfer_vector: &mut FixedVec = &mut [0.0; DIMENTIONALITY];

    let mut data_instance: CSVDataInstance = CSVDataInstance { vel_tan: (0.0), time_1: (0.0), time_2: (0.0), time_3: (0.0), time_4: (0.0), time_5: (0.0) };
    let data_ref: &mut CSVDataInstance = &mut data_instance;
    
    println!("Starting calculations");
    let start = Instant::now();
    for i in (0..=NUM_VELOCITIES).filter(|&i| i != NUM_VELOCITIES/2){
        let mut orbits: i32 = 0;
        let mut time_elapsed: FRes = 0.0;
        let mut theta: FRes = 0.0;
        let tangential_velocity: FRes = MIN_TAN_VELOCITY + V_STEP*{i as FRes};
        println!("index: {i}, v_tan = {tangential_velocity}");
        data_ref.vel_tan = tangential_velocity;

        satellite_ref.velocity = [tangential_velocity, 0.0, 0.0];
        satellite_ref.position = [0.0, ORBITAL_RADIUS, 0.0];
        satellite_ref.acceleration = [0.0, 0.0, 0.0];

        earth_ref.velocity = [0.0, 0.0, 0.0];
        earth_ref.position = [0.0, 0.0, 0.0];
        earth_ref.acceleration = [0.0, 0.0, 0.0]; 

        let mut dump_array: [FRes; ORBITS_TO_COMPLETE as usize] = [0.0; ORBITS_TO_COMPLETE as usize];
        let dump_ref = &mut dump_array;

        //initialising accelerations for improved performance
        set_accelerations(&earth_ref.position, &satellite_ref.position, 
            &mut earth_ref.acceleration, &mut satellite_ref.acceleration, 
            &earth_ref.mass, &satellite_ref.mass, transfer_vector);

        while orbits < ORBITS_TO_COMPLETE{
            let mut old_displacement: FixedVec = satellite_ref.position;
            leap_frog(earth_ref, satellite_ref, transfer_vector);
            theta += get_delta_theta(&mut old_displacement, &mut satellite_ref.position);
            time_elapsed += TIME_STEP;
            if theta.is_nan(){
                println!("Theta is not a number skipping!");
                for i in orbits..ORBITS_TO_COMPLETE{
                    dump_ref[i as usize] = NAN;
                }
                orbits = ORBITS_TO_COMPLETE;
            } else {
                if theta >= 2.0*PI{
                    dump_ref[orbits as usize] = time_elapsed;
                    orbits += 1;
                    theta -= 2.0*PI;
                }
            }
        }
        data_ref.dump_time_array(dump_ref);
        csv_writer.serialize(&data_ref)?;
    }
    let comp_time = start.elapsed();
    println!("time elapsed: {:?}", comp_time);
    csv_writer.flush()?;
    Ok(())
}
