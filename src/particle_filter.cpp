/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
// #include<bits/stdc++.h>

#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100; //Number of particles

	default_random_engine gen;
	double std_x, std_y, std_theta;

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	double weight = 1.0;

	for(int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = weight;
		particle.id = i;
		weights.push_back(weight);
		particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	// Here I have added noise to velocity and yaw_rate based on the noise in x,y and theta.
	// If the results are good then we will keep this, otherwise we will add noise to
	// the final x,y and theta
	// double std_vel = sqrt((pow(std_pos[0],2)) + pow(std_pos[1],2)) / delta_t;
	// double std_yaw_rate = std_pos[2] / delta_t;

	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for(int i = 0; i < num_particles; i++)
	{
		// normal_distribution<double> dist_vel(velocity, std_vel);
		// normal_distribution<double> dist_yaw_rate(yaw_rate, std_yaw_rate);
		// double velocity_sample = dist_vel(gen);
		// double yaw_rate_sample = dist_yaw_rate(gen);
		double x_new;
		double y_new;
		double theta_new;

		if(fabs(yaw_rate < 1e-10))
		{
			x_new = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y_new = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			theta_new = particles[i].theta;
		}

		else
		{
			x_new = particles[i].x + velocity/yaw_rate * 
						(sin(particles[i].theta + yaw_rate * delta_t) - 
						sin(particles[i].theta));
			y_new = particles[i].y + velocity/yaw_rate * 
						(cos(particles[i].theta) - cos(particles[i].theta + 
						yaw_rate * delta_t));
			theta_new = particles[i].theta + yaw_rate * delta_t;
		}
		
		particles[i].x = x_new + dist_x(gen);
		particles[i].y = y_new + dist_y(gen);
		particles[i].theta = theta_new + dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	// cout << "Predicted size: " << predicted.size();
	// cout << "Observation size: " << observations.size();
	for(int i = 0; i < observations.size(); i++)
	{
		double min = 21474654;
		LandmarkObs o = observations[i];
		int predicted_id;
		for(int j = 0; j < predicted.size(); j++)
		{
			// cout << "prediction: " << i << endl;
			LandmarkObs p = predicted[j];

			double diff = dist(o.x, o.y, p.x, p.y);
			if(diff < min)
			{
				min = diff;
				predicted_id = p.id;
				// cout << "observation: " << j << endl;
			}
		}
		observations[i].id = predicted_id;
		// associations.push_back(std::vector<int>{pair[0], pair[1]});
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	for(unsigned int i = 0; i < num_particles; i++)
	{
		// cout << "INSIIIIIIIIIIIIDE" << endl;
		double px = particles[i].x;
		double py = particles[i].y;
		double ptheta = particles[i].theta;
		
		vector<LandmarkObs> landmarks;

		for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			float landmark_x = map_landmarks.landmark_list[j].x_f;
			float landmark_y = map_landmarks.landmark_list[j].y_f;
			float landmark_id = map_landmarks.landmark_list[j].id_i;

			if(fabs(landmark_x - px) <= sensor_range && fabs(landmark_y - py) <= sensor_range)
				landmarks.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
			// cout << landmarks[i].id << endl;
		}

 
	vector<LandmarkObs> global_observation;
	for(unsigned int j = 0; j < observations.size(); j++)
	{

		double t_x = cos(ptheta)*observations[j].x - sin(ptheta)*observations[j].y + px;
		double t_y = sin(ptheta)*observations[j].x + cos(ptheta)*observations[j].y + py;

		global_observation.push_back(LandmarkObs{observations[j].id, t_x, t_y});
	}

	dataAssociation(landmarks, global_observation);

	particles[i].weight = 1.0;
	weights[i] = 1.0;

	for(unsigned int j = 0; j < global_observation.size(); j++)
	{
		double o_x, o_y, pr_x, pr_y;
		int id = global_observation[j].id;
		double s_x = std_landmark[0];
		double s_y = std_landmark[1];

		o_x = global_observation[j].x;
		o_y = global_observation[j].y;

		for(unsigned int k = 0; k < landmarks.size(); k++)
		{
			if(landmarks[k].id == id)
			{
				pr_x = landmarks[k].x;
				pr_y = landmarks[k].y;
			}
		}

		double obs_w = ( 1/(2*M_PI*s_x*s_y)) * 
						exp( -( pow(pr_x-o_x,2)/(2*pow(s_x, 2)) + 
						(pow(pr_y-o_y,2)/(2*pow(s_y, 2))) ) );

		particles[i].weight *= obs_w;
		weights[i] *= obs_w;
	}

	}
}


void ParticleFilter::resample() 
{
	default_random_engine gen;

	uniform_int_distribution<int> uniintdist(0,num_particles-1);
	auto index = uniintdist(gen);

	double max_weight = *max_element(weights.begin(),weights.end());
	uniform_real_distribution<double>uniintdist_(0.0,2*max_weight);

	double beta = 0.0;
	vector<Particle>rmsp_particles;
	for(int i=0;i<num_particles;i++)
	{
		beta+=uniintdist_(gen);
		while(beta>weights[index])
		{
			beta-=weights[index];
			index=(index+1)%num_particles;
		}
		rmsp_particles.push_back(particles[index]);
	}
	particles = rmsp_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
