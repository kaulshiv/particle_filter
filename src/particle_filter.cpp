/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  num_particles = 1000;
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  double sample_x, sample_y, sample_theta;

  for(int i=0; i<num_particles; i++){ 
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
    Particle p = {i, x + sample_x, y + sample_y, theta + sample_theta, 1, {}, {}, {}};
    particles.push_back(p);
  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
  std::default_random_engine gen;
  double sample_x, sample_y, sample_theta, theta0;

  
  for(int i=0; i<num_particles; i++){
    theta0 = particles[i].theta;
    std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
    if(yaw_rate!=0){
      particles[i].x += sample_x + velocity/yaw_rate*(sin(theta0 + yaw_rate*delta_t) - sin(theta0));
      particles[i].y += sample_y + velocity/yaw_rate*(cos(theta0) - cos(theta0 + yaw_rate*delta_t));
      particles[i].theta += sample_theta + yaw_rate*delta_t;
    }
    else{
      particles[i].x += sample_x + velocity*delta_t*cos(theta0);
      particles[i].y += sample_y + velocity*delta_t*sin(theta0);
      particles[i].theta += sample_theta;
    }
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    double min_dist=0;
    int min_id=-1;

    for(int i=0; i<observations.size(); i++){
      double min_dist=0;
      int min_id=-1;
      // x = observations[i].x;
      // y = observations[i].y;



      for(int j=0; j<predicted.size(); j++){
        printf("test");
      }
    }



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double x, y, xt, yt, transformed_x, transformed_y, theta;

  for(int i=0; i<num_particles; i++){
    vector<LandmarkObs> predicted = {};
    for(int j=0; j<map_landmarks.landmark_list.size(); j++){
      if(dist(particles[i].x, particles[i].y, (double)map_landmarks.landmark_list[j].x_f, 
         (double)map_landmarks.landmark_list[j].y_f) <= sensor_range){
          
        x = map_landmarks.landmark_list[j].x_f;
        y = map_landmarks.landmark_list[j].y_f;
        xt = particles[i].x;
        yt = particles[i].y;
        theta = particles[i].theta;

        transformed_x = x*cos(theta) - y*sin(theta) + xt;
        transformed_y = x*sin(theta) + y*cos(theta) + yt;
        LandmarkObs obs = {map_landmarks.landmark_list[j].id_i, transformed_x, transformed_y};
        predicted.push_back(obs);
      }
    }

    dataAssociation(predicted, observations);



   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
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

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}