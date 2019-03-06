/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 * Modified on: Mar 06, 2019
 * Author: Simon Wisseler
 */

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta) {
    /**
     * Set number of particles and initialize each particle to first position
     * based on (noisy) GPS data
     */
    num_particles = 25;
    
    for (int i = 0; i < num_particles; i++) {
        Particle p;
        p.id = i;
        p.x = x;
        p.x += gaussian_noise_x(gen);
        p.y = y;
        p.y += gaussian_noise_y(gen);
        p.theta = theta;
        p.theta += gaussian_noise_theta(gen);
        p.weight = 1.0;
        particles.push_back(p);
    }
    is_initialized = true;
}


void ParticleFilter::predictGlobalCoord(double delta_t, double v, double w) {
    /**
     * Update particles using motion model and gaussian noise
     * NOTE: Not applying sample_motion_model_velocity from Probabilistic Robotics - S. Thrun, W. Burgard, D. Fox (2006), p. 124, as
     * no standard deviation for velocity and yaw rate provided in main.cpp
     */
    for (int i = 0; i < num_particles; i++) {
        
        if (fabs(w) < 0.00001) {
            particles[i].x += v * delta_t * cos(particles[i].theta);
            particles[i].y += v * delta_t * sin(particles[i].theta);
        }
        else {
            particles[i].x += v / w * (sin(particles[i].theta + w*delta_t) - sin(particles[i].theta));
            particles[i].y += v / w * (cos(particles[i].theta) - cos(particles[i].theta + w*delta_t));
            particles[i].theta += w * delta_t;
        }
        particles[i].x += gaussian_noise_x(gen);
        particles[i].y += gaussian_noise_y(gen);
        particles[i].theta += gaussian_noise_theta(gen);
    }
}

vector<LandmarkObs> ParticleFilter::predictMeasurements(Particle p, const Map &map_lms){
    /**
     * For a given particle, predict measurements from particle to landmarks within sensor range
     */
    vector<LandmarkObs> predictions;
    
    for (unsigned int j = 0; j < map_lms.landmark_list.size(); j++) {
        double lm_x = map_lms.landmark_list[j].x_f;
        double lm_y = map_lms.landmark_list[j].y_f;
        int lm_id = map_lms.landmark_list[j].id_i;
        
        if (dist(lm_x, lm_y, p.x, p.y) <= sensor_range) {
            predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
        }
    }
    return predictions;
};

vector<LandmarkObs> ParticleFilter::transformObsFromVehicleCoord2GlobalCoord(Particle p, const vector<LandmarkObs> &obs){
    /**
     * Transform sensor measurements from vehicle's local coordinate system to global coordinate system
     */
    vector<LandmarkObs> transformed_obs;
    
    for (unsigned int i = 0; i < obs.size(); i++) {
        double transformed_x = cos(p.theta)*obs[i].x - sin(p.theta)*obs[i].y + p.x;
        double transformed_y = sin(p.theta)*obs[i].x + cos(p.theta)*obs[i].y + p.y;
        transformed_obs.push_back(LandmarkObs{ obs[i].id, transformed_x, transformed_y });
    }
    return transformed_obs;
};

void ParticleFilter::associateData(const vector<LandmarkObs> &pred, vector<LandmarkObs> &obs) {
    /**
     * For a given particle, compare predicted measurements with observed measurements from sensors and
     * associate observed measurements with closest landmark
     */
    for (unsigned int i = 0; i < obs.size(); i++) {
        LandmarkObs o = obs[i];
        // Initialize variables
        double min_dist = numeric_limits<double>::max();
        int map_id = -1;
        // Loop through predicted measurements to determine closest landmark
        
        for (unsigned int j = 0; j < pred.size(); j++) {
            LandmarkObs p = pred[j];
            double cur_dist = dist(o.x, o.y, p.x, p.y); //dist() defined in helper_functions.h
            
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                map_id = p.id;
            }
        }
        obs[i].id = map_id;
    }
}

void ParticleFilter::updateWeights(double std_landmark[2], const vector<LandmarkObs> &obs, const Map &map_lms) {
    /**
     * Update weight of each particle
     */
    for (unsigned int i = 0; i < num_particles; i++) {
        const vector<LandmarkObs> pred = predictMeasurements(particles[i], map_lms);
        vector<LandmarkObs> transformed_obs = transformObsFromVehicleCoord2GlobalCoord(particles[i], obs);
        associateData(pred, transformed_obs);
        particles[i].weight = 1.0; // (re)initialize weight
        
        for (unsigned int j = 0; j < transformed_obs.size(); j++) {
            double obs_x, obs_y, pred_x, pred_y;
            obs_x = transformed_obs[j].x;
            obs_y = transformed_obs[j].y;
            int associated_pred = transformed_obs[j].id;
            
            for (unsigned int k = 0; k < pred.size(); k++) {
                if (pred[k].id == associated_pred) {
                    pred_x = pred[k].x;
                    pred_y = pred[k].y;
                }
            }
            double s_x = std_landmark[0];
            double s_y = std_landmark[1];
            double obs_w = (1/(2*M_PI*s_x*s_y)) * exp( -(pow(pred_x-obs_x,2)/(2*pow(s_x, 2)) + (pow(pred_y-obs_y,2)/(2*pow(s_y, 2)))));
            particles[i].weight *= obs_w;
        }
    }
}

void ParticleFilter::resample() {
    /**
     * Resample particles with replacement with probability proportional to their respective weights
     * using resampling wheel (see https://www.youtube.com/watch?v=wNQVo6uOgYA)
     */
    vector<Particle> new_particles;
    vector<double> weights;
    
    for (unsigned int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    uniform_int_distribution<int> uni_int_dist(0, num_particles-1);
    auto index = uni_int_dist(gen); // random starting index
    double max_weight = *max_element(weights.begin(), weights.end());
    uniform_real_distribution<double> uni_real_dist(0.0, max_weight);
    double beta = 0.0;
    
    for (unsigned int i = 0; i < num_particles; i++) {
        beta += uni_real_dist(gen) * 2.0;
        
        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        new_particles.push_back(particles[index]);
    }
    particles = new_particles;
}


void ParticleFilter::SetAssociations(Particle& particle, const vector<int>& associations, const vector<double>& sense_x, const vector<double>& sense_y) {
    /**
     * particle: the particle to which assign each listed association,
     * and association's (x,y) world coordinates mapping
     * associations: The landmark id that goes along with each listed association
     * sense_x: the associations x mapping already converted to world coordinates
     * sense_y: the associations y mapping already converted to world coordinates
     */
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}


string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseCoord(Particle best, string coord) {
    vector<double> v;
    
    if (coord == "X") {
        v = best.sense_x;
    } else {
        v = best.sense_y;
    }
    
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
