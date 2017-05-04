#ifndef SWARM_INDIVIDUAL__h
#define SWARM_INDIVIDUAL__h

#include <vector>
using namespace std;

class SwarmIndividual {
public:
	SwarmIndividual();
	virtual ~SwarmIndividual();

	void SetParams(vector<double> P) { Params = P; };
	vector<double> GetParams() { return Params; };

	void SetVelocity(vector<double> V) { Velocity = V; };
	vector<double> GetVelocity() { return Velocity; };

	void SetLocalBestParams(vector<double> P) { LocalBestParams = P; };
	vector<double> GetLocalBestParams() { return LocalBestParams; };

	void SetLocalBest(double b) { LocalBestValue = b; };
	double GetLocalBest() { return LocalBestValue; };

private:
	vector<double> Params;
	vector<double> Velocity;
	vector<double> LocalBestParams;
	double LocalBestValue;

};

#endif
