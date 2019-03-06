/**
 * particle_filter.h
 * 2D particle filter class.
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 * Modified on: Mar 06, 2019
 * Author: Simon Wisseler
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using namespace std;

struct Particle {
    int id;
    double x;
    double y;
    double theta;
    double weight;
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
};


class ParticleFilter {
public:
    /**
     * Constructor
     * @param arr Array of standard deviations (x, y, theta) required for generation of gaussian noise via normal distributions
     * @param s_range Sensor_range [m]
     */
    ParticleFilter(double (&arr)[3], double sensor_range);
    
    //Destructor
    ~ParticleFilter();
    
    /**
     * init Initialize each particle to first position based on (noisy) GPS data
     * @param x Initial x position [m]
     * @param y Initial y position [m]
     * @param theta Initial orientation [rad]
     */
    void init(double x, double y, double theta);
    
    /**
     * predictGlobalCoord Update particles using motion model and gaussian noise
     * @param delta_t Time between time step t and t+1 in measurements [s]
     * @param v Velocity of car from t to t+1 [m/s]
     * @param w Yaw rate of car from t to t+1 [rad/s]
     */
    void predictGlobalCoord(double delta_t, double v, double w);
    
    /**
     * predictMeasurements For a given particle, predict measurements from particle to landmarks within sensor range
     * @param p Particle of interest
     * @param map_lms Map class containing all map landmarks
     * @return predictions Vector of predicted distances between particle of interest and all landmarks within sensor range
     */
    vector<LandmarkObs> predictMeasurements(Particle p, const Map &map_lms);
    
    /**
     * transformObsFromVehicleCoord2GlobalCoord Transform sensor measurements from vehicle's local coordinate system to global coordinate system
     * @param p Particle of interest
     * @param observations Vector of landmark observations in vehicle's local coordinate system
     * @return transformed_obs Vector of landmark observations in global coordinate system
     */
    vector<LandmarkObs> transformObsFromVehicleCoord2GlobalCoord(Particle p, const vector<LandmarkObs> &observations);
    
    /**
     * associateData For a given particle, compare predicted measurements with observed measurements from sensors and
     * associate observed measurements with closest landmark
     * @param pred Vector of predicted distances between particle of interest and all landmarks within sensor range
     * @param obs observations Vector of landmark observations in global coordinate system
     */
    void associateData(const vector<LandmarkObs> &pred, vector<LandmarkObs> &obs);
    
    /**
     * updateWeights Update weight of each particle
     * @param std_landmark[] Standard deviaton of landmark measurements [x [m], y [m]]
     * @param observations Vector of landmark observations in global coordinate system
     * @param map Map class containing all map landmarks
     */
    void updateWeights(double std_landmark[2],
                       const vector<LandmarkObs> &observations,
                       const Map &map_landmarks);
    
    /**
     * resample Resample particles with replacement with probability proportional to their respective weights
     * using resampling wheel (see https://www.youtube.com/watch?v=wNQVo6uOgYA)
     */
    void resample();
    
    /**
     * Set a particles list of associations, along with the associations'
     *   calculated world x,y coordinates
     * This can be a very useful debugging tool to make sure transformations
     *   are correct and assocations correctly connected
     */
    void SetAssociations(Particle& particle, const vector<int>& associations,
                         const vector<double>& sense_x,
                         const vector<double>& sense_y);
    
    /**
     * initialized Returns whether particle filter is initialized yet or not
     */
    const bool initialized() const {
        return is_initialized;
    }
    
    /**
     * Used for obtaining debugging information related to particles
     */
    string getAssociations(Particle best);
    string getSenseCoord(Particle best, string coord);
    
    // Set of current particles
    vector<Particle> particles;
    
private:
    // Number of particles to draw
    int num_particles;
    
    // Flag, if filter is initialized
    bool is_initialized;
    
    // Vector of weights of all particles
    vector<double> weights;
    
    // Array of standard deviations (x, y, theta) required for generation of gaussian noise via normal distributions
    double (&std_pos)[3];
    
    // Normal distributions from which to sample noise for x, y, theta coordinate
    normal_distribution<double> gaussian_noise_x;
    normal_distribution<double> gaussian_noise_y;
    normal_distribution<double> gaussian_noise_theta;
    
    // Maximum sensor range
    double sensor_range;
};

#endif  // PARTICLE_FILTER_H_
