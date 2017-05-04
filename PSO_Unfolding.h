#ifndef PSO_UNFOLDING__h
#define PSO_UNFOLDING__h

#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream> 
#include <string>
#include <random>

using namespace std; 

const int Nfoils = 15;
const int NEBins = 100;
vector<double> Resp[Nfoils];
double act_measured[Nfoils] = { ADD MEASUREMENT }; // !! removed measurement before sharing code

vector<double> ReadResponseFunction(const char* FileName, int foilNo, const int Nbins) {
	vector<double> v;
	ifstream in(FileName);
	//ignore header
	in.ignore(10000, '\n'); 
	string s; int nFoil; double x;
	while (1) {
		if (!in.good())
			break;

		in >> s >> nFoil;
		if (nFoil != (foilNo+1) ) {
			in.ignore(100000000, '\n');
		}
		else {
			for (int i = 0; i < Nbins; i++) {
				in >> x;
				v.push_back(x);
			}
			return v;
		}
	}
	return v;
}


/*******************************************************************************/
/* Function updating position and velocity values:
v_new = w0 v_old + w1 r1 (localBestPos - p_old) + w2 r2 (globalBest - p_old)
with w0, w1, w2: weights
r1, r2: random numbers in [0,1]
p_new = p_old + v_new                                                       */
/*******************************************************************************/
static void UpdateParameters(SwarmIndividual &ThisBee, vector<double> GlobBest, const double w0_min, const double w0_max, const double w1_min, const double w1_max, const double w2_min, const double w2_max, const int nSteps, const int current) {
	double r1, r2, v_old, p_old, w0_i, w1_i, w2_i, p_new;

	/* time dependence */
	w0_i = (w0_max - w0_min)*(nSteps - current) / nSteps + w0_min;
	w1_i = (w1_min - w1_max)*current / nSteps + w1_max;
	w2_i = (w2_max - w2_min)*current / nSteps + w2_min;

	int nPar = (int)ThisBee.GetParams().size();
	vector<double> newParams;
	vector<double> newVelocity;
	for (int h = 0; h < nPar; h++) {
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<> dis(0, 1);
		r1 = dis(gen);
		r2 = dis(gen);
		if (r1 == r2)
			cout << "ERROR in random numbers: update params" << endl;

		v_old = ThisBee.GetVelocity().at(h);
		p_old = ThisBee.GetParams().at(h);

		newVelocity.push_back(w0_i*v_old + w1_i*r1*(ThisBee.GetLocalBestParams().at(h) - p_old) + w2_i*r2*(GlobBest.at(h) - p_old));

		/* no boarders yet*/
		/*
		v_max = 0.5*fabs(Pmax[h] - Pmin[h]);
		if (fabs(ThisBee[BeeNo].Velocity[h]) > v_max)
			ThisBee[BeeNo].Velocity[h] = ((ThisBee[BeeNo].Velocity[h] > 0) ? v_max : -v_max);
		*/
		p_new = p_old + newVelocity.at(h);
		if (p_new < 0) {
			p_new = 0;
		}
		// set max. diff. between bins to: factor 50 
		if (h > 0 && p_new > newParams.at(h - 1) * 50 ) {
			p_new = newParams.at(h - 1) * 50;
		}
		else if (h > 0 && p_new < newParams.at(h - 1) / 50) {
			p_new < newParams.at(h - 1) / 50;
		}

		newParams.push_back(p_new);

	}
	ThisBee.SetParams(newParams);
	ThisBee.SetVelocity(newVelocity);
}


/********************************/
// update best FoMs 
/*******************************/
double UpdateFoM(vector<double> Params, int FoilCount) {
	double FoM = 0;
	double act = 0;
	double weight = 0;
	// TEST: give higher weight to first foils (1/2)------------------
	double weightSum = 0;
	for (int i = 0; i < FoilCount; i++) {
		weightSum += Nfoils - i;
	}
	//-------------------------------------------------------------
	const int nParams = (int)Params.size();
	for (int f = 0; f < Nfoils; f++) {
		weight = (FoilCount > f) ? (double)(1.0 / FoilCount) : 0;
		if (!weight)
			break;
		// TEST: give higher weight to first foils (2/2)------------------
		weight = (15 - f) / weightSum;
		//--------------------------------------------------------------
		act = 0;
		for (int e = 0; e < nParams; e++) {
			act += Params.at(e) * Resp[f].at(e);
		}
		FoM += weight * (act - act_measured[f])*(act - act_measured[f]); 
	}
	return sqrt(FoM);
}


#endif
