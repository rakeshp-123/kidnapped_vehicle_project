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
#include <float.h>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;
	default_random_engine gen;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for(int i = 0; i < num_particles;  i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(p.weight);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	double x0 = 0.0;
	double y0 = 0.0;
	double theta0 = 0.0;

	default_random_engine gen;

	for(int i = 0; i < num_particles;  i++)
	{
		Particle p;
		x0 = particles[i].x;
		y0 = particles[i].y;
		theta0 = particles[i].theta;
		
		//Checking for near zero yaw rate and applying appropriate motion equation based on that.
		if(fabs(yaw_rate) >= 0.0001)
		{
			p.x = x0 + ((velocity/yaw_rate) * (sin(theta0 + (yaw_rate * delta_t)) - sin(theta0)));
			p.y = y0 + ((velocity/yaw_rate) * (cos(theta0) - cos(theta0 + (yaw_rate * delta_t))));
			p.theta = theta0 + (yaw_rate * delta_t);
		}

		else
		{
			p.x = x0 + velocity * cos(theta0) * delta_t;
			p.y = y0 + velocity * sin(theta0) * delta_t;
			p.theta = theta0;
		}

		normal_distribution<double> dist_x(p.x, std_pos[0]);
		normal_distribution<double> dist_y(p.y, std_pos[1]);
		normal_distribution<double> dist_theta(p.theta, std_pos[2]);
		
		//initializing 10 partiacles at random within provided normalized distribution range.
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double minDist;
	double obs_x;
	double obs_y;
        double distCurr;
	int id;

	for(int i = 0; i < observations.size(); i++)
	{
		minDist = DBL_MAX;
		obs_x = observations[i].x;
		obs_y = observations[i].y;
		id = -1;

		// find the minimum distance of observations to the landmarks.
		for(int j = 0; j < predicted.size(); j++)
		{
			distCurr = dist(obs_x, obs_y, predicted[j].x, predicted[j].y);

			if(distCurr < minDist)
			{
				minDist = distCurr;
				id = predicted[j].id;
			}
		}
		observations[i].id = id;
	}
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html 
	double ym = 0.0;
	double xm = 0.0;
	double thetam = 0.0;

	
	double total_weight = 0.0;

	for(int i = 0; i < num_particles ; i++)
	{
		double xp = particles[i].x;
		double yp = particles[i].y;
		double thetap = particles[i].theta;
		std::vector<LandmarkObs> transformed_obs;

		//transforming the observations from vehicle co-ordinate to map co-ordinate
		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs trnsfrm_ob;
			trnsfrm_ob.id = j;
			trnsfrm_ob.x = xp + (cos(thetap) * observations[j].x) - (sin(thetap) * observations[j].y);
			trnsfrm_ob.y = yp + (sin(thetap) * observations[j].x) + (cos(thetap) * observations[j].y);
			transformed_obs.push_back(trnsfrm_ob);
		}

		//finding the landmarks within acceptable range of a particle
		std::vector<LandmarkObs> particle_to_landmark_pred;
		
		for(int j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
			
			if((fabs(xp - landmark.x_f) <= sensor_range) && (fabs(xp - landmark.x_f) <= sensor_range))
			{
				LandmarkObs pred_lmrk;
				pred_lmrk.id = landmark.id_i;
				pred_lmrk.x = landmark.x_f;
				pred_lmrk.y = landmark.y_f;
				particle_to_landmark_pred.push_back(pred_lmrk);
			}			
		}

		//find the landmark with minimum distance to observations.
		dataAssociation(particle_to_landmark_pred, transformed_obs);


		// assign weight to particles based on the probability using multi-variate probability distribution.
		double sigmaX = std_landmark[0];
		double sigmaY = std_landmark[1];
		double sigmaX2 = pow(sigmaX, 2);
		double sigmaY2 = pow(sigmaY, 2);

		double nor = (1.0/(2.0 * M_PI * sigmaX * sigmaY));

		particles[i].weight = 1.0;

		for(int j=0; j< transformed_obs.size(); j++)
		{
			double prob = 1.0;
			double obs_x = transformed_obs[j].x;
			double obs_y = transformed_obs[j].y;
			int obs_id = transformed_obs[j].id;


			for(int h=0; h< particle_to_landmark_pred.size(); h++)
			{	
				double pred_x = particle_to_landmark_pred[h].x;
				double pred_y = particle_to_landmark_pred[h].y;
				//cout<<"in here/n";
				if(particle_to_landmark_pred[h].id == transformed_obs[j].id)
				{
					prob = nor * exp(-1.0 * ((pow((obs_x - pred_x), 2)/(2.0 * sigmaX2)) + (pow((obs_y - pred_y), 2)/(2.0 * sigmaY2))));
					particles[i].weight *= prob;
				}
			}
		}
		total_weight += particles[i].weight;
	}

	
	//Normalize weights
	for(int k = 0; k < particles.size(); k++)
	{
		particles[k].weight = particles[k].weight/total_weight;
		weights[k] = particles[k].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	

	std::vector<Particle> sampled_particles;

	default_random_engine gen;

	discrete_distribution<int> idxDist(0, num_particles - 1);

	int idx = idxDist(gen);
	double beta = 0.0;
	double max_weight = *max_element(weights.begin(), weights.end());

	//resample particles
	for(int i = 0; i < particles.size(); i++)
	{
		uniform_real_distribution<double> idxDistBeta(0.0, max_weight * 2.0);
		beta += idxDistBeta(gen);
		while(beta  > weights[idx])
		{
			beta -= weights[idx];
			idx = (idx + 1) % num_particles;
		}
		sampled_particles.push_back(particles[idx]);
	}
	particles = sampled_particles;
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
	
    return particle;
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
