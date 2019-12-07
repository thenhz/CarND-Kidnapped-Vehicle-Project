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

void ParticleFilter::init(double x, double y, double theta, double std[])
{
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50; // TODO: Set the number of particles
  weights.resize(num_particles);
  particles.resize(num_particles);

  std::default_random_engine gen;
  std::normal_distribution<double> 
      dist_x(x, std[0]),
      dist_y(y, std[1]),
      dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i)
  {

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0.0, std_pos[0]),
      dist_y(0.0, std_pos[1]),
      dist_theta(0.0, std_pos[2]);

  for (auto &particle : particles)
  {
    if (fabs(yaw_rate) > 1e-5)
    {
      const auto theta_pred = particle.theta + (yaw_rate * delta_t);
      particle.x += (velocity / yaw_rate) * (+sin(theta_pred) - sin(particle.theta));
      particle.y += (velocity / yaw_rate) * (-cos(theta_pred) + cos(particle.theta));
      particle.theta += yaw_rate * delta_t;
    }
    else
    {
      particle.x += velocity * cos(particle.theta) * delta_t;
      particle.y += velocity * sin(particle.theta) * delta_t;
    }

    // add sensor noise
    particle.x += dist_x(gen);
    particle.y += dist_y(gen);
    particle.theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for(int i = 0 ; i<observations.size();i++){
    vector<double> distances;
    
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
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
  // penalty when observation is out of range

    // denominator term
  const double gauss_x_den = 2 * pow(std_landmark[0], 2);
  const double gauss_y_den = 2 * pow(std_landmark[1], 2);
  const double gauss_den   = 2 * M_PI * std_landmark[0] * std_landmark[1];
  
  const double OUT_OF_RANGE_PENALTY = 9999999.999;

  // define associated landmark parameters
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;

  // updated weght for each particle
  double updated_weight = 1.0;

  // iterate particles
  for (auto &particle : particles)
  {

    // reset values
    updated_weight = 1.0;
    associations = {};
    sense_x = {};
    sense_y = {};

    // iterate observations for each particle
    for (auto const &obs_vcs : observations)
    {

      /**************************************************************
       * STEP - 1:
       * Transformation observation
       * From vehicle coordinate system (VCS) to map coordinate system (MCS)
       **************************************************************/

      double obs_mcs_x = obs_vcs.x * cos(particle.theta) - obs_vcs.y * sin(particle.theta) + particle.x;
      double obs_mcs_y = obs_vcs.x * sin(particle.theta) + obs_vcs.y * cos(particle.theta) + particle.y;

      /**************************************************************
       * STEP - 2:
       * Nearest Landmark from Observation
       **************************************************************/

      // find distances of an MCS observation to all map landmarks
      vector<double> distances;
      for (auto const &l : map_landmarks.landmark_list)
      {
        const double dx = pow(particle.x - l.x_f, 2);
        const double dy = pow(particle.y - l.y_f, 2);
        const double distance_from_particle = sqrt(dx + dy);
        double distance = 0;

        // distance between lidar and landmark
        if (distance_from_particle <= sensor_range)
        {
          const double obs_dx = pow(obs_mcs_x - l.x_f, 2);
          const double obs_dy = pow(obs_mcs_y - l.y_f, 2);
          distance = sqrt(obs_dx + obs_dy);
        }

        // penalize observations which are out of sensor range
        else
        {
          distance = OUT_OF_RANGE_PENALTY;
        }

        // append distance
        distances.push_back(distance);
      }

      // associated landmark for MCS observation with minimum distance
      const auto index = indexOfSmallestElement(distances);
      const auto associated_landmark = map_landmarks.landmark_list[index];

      /**************************************************************
       * STEP - 3:
       * Find normalized probability using multi-variate Gaussian distribution
       **************************************************************/

      // argument of exponential term
      double exp_arg = 0.0;
      exp_arg += pow(obs_mcs_x - associated_landmark.x_f, 2) / gauss_x_den;
      exp_arg += pow(obs_mcs_y - associated_landmark.y_f, 2) / gauss_y_den;

      // update weights with normalization of all observations
      updated_weight *= exp(-exp_arg) / gauss_den;

      //updated_weight = multiv_prob(std_landmark[0], std_landmark[1], obs_mcs_x, obs_mcs_y, associated_landmark.x_f, associated_landmark.y_f);

      // append associated landmark
      associations.push_back(associated_landmark.id_i);
      sense_x.push_back(associated_landmark.x_f);
      sense_y.push_back(associated_landmark.y_f);
    }

    /**************************************************************
     * STEP - 4:
     * Update Particle Weight
     **************************************************************/

    SetAssociations(particle, associations, sense_x, sense_y);

    // update particle weight
    particle.weight = updated_weight;

    // index of particle (range-based for)
    const auto index = &particle - &particles[0];

    // update weight
    weights[index] = updated_weight;
  }
}

void ParticleFilter::resample()
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  vector<Particle> resampled_particles(num_particles);
  for (auto &particle : resampled_particles)
  {
    std::discrete_distribution<int> sample_index(weights.begin(), weights.end());
    particle = particles[sample_index(gen)];
  }

  particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}