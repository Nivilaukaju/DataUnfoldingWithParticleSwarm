#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include "SwarmIndividual.h"
#include <vector>
#include <time.h>
#include <random>
#include <algorithm> 
#include <functional>
#include <string>
#include "PSO_Unfolding.h"

using namespace std;

int main(int argc, char* argv[]) {


	/* swarm parameters: user input */
	if (argc != 4) {
		cout << "ERROR: need input arguments:" << endl
			<< " (1) number of swarm individuals" << endl
			<< " (2) number of parameters (e-bins)" << endl
			<< " (3) number of iterations" << endl
			<< endl << "- press any key to terminate and try again" << endl;
		cin.get();
		return -1;
	}
	int Nbees_tmp = (int)atof(argv[1]);
	int Npar = (int)atof(argv[2]);
	int Nit = (int)atof(argv[3]);
	//overwrite values for testing:
	const int Nbees = 10000;
	//Nit = 150;
	//Npar = 100;

	const int NitBreak = Nit; //not used right now
	const int MaxItWithoutImprovement = 20; //a) break, b) move swarm indiv. --> right now: a)

	/* input parameters with default values */
	double wMin_inertia = 0.4,      // min. inertia weight (default: 0.4)
		wMax_inertia = 0.9,         // max. inertia weight (default: 0.9)
		wMin_local = 0.5,           // min. weight individual best (default: 0.5)
		wMax_local = 2.5,           // max. weight individual best (default: 2.5)
		wMin_global = 0.5,          // min. weight global best (default: 0.5)
		wMax_global = 2.5;          // max. weight global best (default: 2.5)


	/* global best values*/
	double GlobalBestFoM = 1e12; //minimize this
	vector<double> GlobalBestParams;

	/* write log-file*/
	FILE *fOutfile;
	const char* sOutfileName = "MyOutfile.txt";
	fOutfile = fopen(sOutfileName, "wt");
	if (!fOutfile) {
		fprintf(stderr, "opening %s: %s\n", sOutfileName, strerror(errno)); 
		return 1;
	}
	fprintf(fOutfile, "Particle Swarm Optimization logfile: \n\n");
	fprintf(fOutfile, "User input parameters : \n");
	fprintf(fOutfile, " Number of swarm individuals              : %d \n", Nbees);
	fprintf(fOutfile, " Number of parameters (output energy bins): %d \n", Npar);
	fprintf(fOutfile, " Nominal number of iterations             : %d \n", Nit);
	fprintf(fOutfile, " Actual number of iterations              : %d \n", NitBreak);
	fprintf(fOutfile, " Maximum number of iterations without improvement: %d \n", MaxItWithoutImprovement);
	fprintf(fOutfile, "\n Default parameters : \n");
	fprintf(fOutfile, " inertia weights: varying between %f and %f \n", wMin_inertia, wMax_inertia);
	fprintf(fOutfile, " local weights  : varying between %f and %f \n", wMin_local, wMax_local);
	fprintf(fOutfile, " global weights : varying between %f and %f \n", wMin_global, wMax_global);


	//get ResponseFunctions:
	vector<double> Resp0[Nfoils];
	for (int f = 0; f < Nfoils; f++) {
		Resp0[f] = ReadResponseFunction("Path/ResponseFunctions.fmt", f, NEBins);

		if (NEBins != Npar) {
			if (NEBins%Npar != 0) {
				cout << "ERROR: Npar must be a divider of 100" << endl;
			}
			const int nRebin = NEBins / Npar;
			if (nRebin == 0) {
				cout << "ERROR: nRebin=0" << endl;
				break;
			}
			for (int e = 0; e < NEBins - (nRebin - 1); e++) {
				double tmp = 0;
				for (int i = 0; i < nRebin; i++) {
					tmp += Resp0[f].at(e + i);
				}
				tmp /= nRebin;
				e += nRebin - 1;
				Resp[f].push_back(tmp);
			}
		}
		else {
			Resp[f] = Resp0[f];
		}
		if ((int)Resp[f].size() != Npar)
			cout << "ERROR: some bug in rebin: size Resp=" << (int)Resp[f].size() << ", npar=" << Npar << endl;
	}//foils

	//================================================
	SwarmIndividual *OneBee = new SwarmIndividual[Nbees]; 
	/* initialize one bee with Odin spectrum shape: */
	double ParOdinNormalized[50]={ 3.10E-04,3.10E-04,3.10E-04, 5.48E-04, 9.73E-04, 1.72E-03, 3.05E-03, 5.40E-03, 2.01E-02, 5.09E-02, 1.21E-01, 2.28E-01, 2.87E-01, 1.85E-01, 5.27E-02, 4.73E-03, 2.00E-03, 1.72E-03, 1.53E-03, 1.39E-03, 1.05E-03, 1.28E-03, 1.08E-03, 1.06E-03, 1.05E-03, 8.99E-04, 8.02E-04, 9.94E-04, 7.90E-04, 8.95E-04, 7.72E-04, 7.15E-04, 6.78E-04, 6.41E-04, 5.88E-04, 5.42E-04, 4.76E-04, 4.36E-04, 3.92E-04, 3.64E-04, 4.59E-04, 1.00E-03, 2.55E-03, 7.18E-03, 5.60E-03,5.60E-03,5.60E-03,5.60E-03,5.60E-03,5.60E-03 };
	for (int i = 0; i < Nbees; i++) {
		vector<double> p0, v0;
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<> dis(-6, 6);
		double k = pow(10, dis(gen));
		for (int j = 0; j < Npar; j++) {
		if (i == 0) {
				k = 3.88537e6; //first bee at simulated result (scale) 
			}
			double ParSim2 = ParOdinNormalized[(int)(j/2)]; //for Npar=100
			if(Npar==50) 
				ParSim2 = ParOdinNormalized[j];
			p0.push_back(k*ParSim2);
			v0.push_back(0.1*k*ParSim2);
		}
		OneBee[i].SetParams(p0);
		OneBee[i].SetLocalBestParams(p0);
		OneBee[i].SetVelocity(v0);

		//calc figure of merit (FoM): 
		// sum of activation difference (calc-meas) in each foil
		double FoM = UpdateFoM(OneBee[i].GetParams(), 1); 
		OneBee[i].SetLocalBest(FoM);
		if (FoM < GlobalBestFoM) {
			GlobalBestFoM = FoM;
			GlobalBestParams = OneBee[i].GetLocalBestParams();
		}
	}

	fprintf(fOutfile, "swarm initialized: best staring value \n");
	fprintf(fOutfile, " global best FoM: %f \n",GlobalBestFoM);
	//DEBUG
	fprintf(fOutfile, " at parameters: \n");
	for (int p = 0; p < Npar; p++) {
		fprintf(fOutfile, "%10.3f   ", GlobalBestParams.at(p));
	}

	//================================================
	// Optimization loop
	//================================================
	double FoM = 0, act=0;
	int NItWoImprovement = 0;
	int FoM_weightChange = NitBreak / Nfoils;
	int FoM_weightCount = 1; 
	for (int it = 0; it < NitBreak; it++) {
		cout << "... working on iteration " << it << ", current best FoM: " << GlobalBestFoM << endl;
		bool FoundBetterFoM = false;
		if (it > FoM_weightChange*FoM_weightCount && FoM_weightCount<Nfoils ) {
			FoM_weightCount++;
			// re-calc current best FoM including more foils:
			GlobalBestFoM=UpdateFoM(GlobalBestParams, FoM_weightCount);

			for (int b = 0; b < Nbees; b++) {
				OneBee[b].SetLocalBest(UpdateFoM(OneBee[b].GetLocalBestParams(), FoM_weightCount));
			}
		}
		for (int b = 0; b < Nbees; b++) {
			UpdateParameters(OneBee[b], GlobalBestParams, wMin_inertia, wMax_inertia, wMin_local, wMax_local, wMin_global, wMax_global, Nit, it);
			FoM = UpdateFoM(OneBee[b].GetParams(), FoM_weightCount);
			if (FoM < OneBee[b].GetLocalBest()) {
				OneBee[b].SetLocalBest(FoM);
				OneBee[b].SetLocalBestParams(OneBee[b].GetParams());
				if (FoM < GlobalBestFoM) {
					GlobalBestFoM = FoM;
					GlobalBestParams = OneBee[b].GetParams();
					FoundBetterFoM = true;
				}
			}
		}//bees
		if (FoundBetterFoM)
			NItWoImprovement = 0;
		else if (!FoundBetterFoM) {
			NItWoImprovement++;
			// what to do if no improvement for X iterations:
			if (NItWoImprovement > MaxItWithoutImprovement) {
				//a) break;
				// break;
				
				//b) move some individuals
				/*for (int b = 0; b < Nbees; b++) {
					random_device rd;
					mt19937 gen(rd());
					uniform_real_distribution<> dis(-1, 1);
					double MyRnd = dis(gen);
					if (fabs(MyRnd) > 0.5) {
						// re-ini:
						double SumFlux = 0;
						for (int j = 0; j < Npar; j++) {
							SumFlux += GlobalBestParams.at(j);
						}
						double MeanFlux = SumFlux / Npar;
						vector<double> newp0, newv0;
						random_device rd2;
						mt19937 gen2(rd2());
						uniform_real_distribution<> dis2(MeanFlux*0.8, MeanFlux*1.2);
						double k2 =  dis2(gen2);
						for (int j = 0; j < Npar; j++) {
							newp0.push_back(k2);
							newv0.push_back(0.1*k2);
						}
						OneBee[b].SetParams(newp0);
						OneBee[b].SetVelocity(newv0);
						//---------------------------------------
					}
				} */

				// a') change FoM to final before break next time
				if (FoM_weightCount < Nfoils) {
					FoM_weightCount = Nfoils;
					NItWoImprovement = 0;
					GlobalBestFoM = UpdateFoM(GlobalBestParams, FoM_weightCount);
					cout << "set counting foils to all; updated FoM: " << GlobalBestFoM << endl;
					for (int b = 0; b < Nbees; b++) {
						OneBee[b].SetLocalBest(UpdateFoM(OneBee[b].GetLocalBestParams(), FoM_weightCount));
					}
					continue;
				}
				else
					break;
				NItWoImprovement = 0;
			}
		}


		fprintf(fOutfile, "\n \n (%d) \n",it);
		fprintf(fOutfile, " global best FoM: %f \n", GlobalBestFoM);
		fprintf(fOutfile, " at parameters: \n");
		for (int p = 0; p < Npar; p++) {
			fprintf(fOutfile, "%10.3f   ", GlobalBestParams.at(p));
		}

	}//iterations
	fclose(fOutfile);
	//=============================
	// cleanup
	cout << "final global best FoM = " << GlobalBestFoM << endl<<"  at params: "<<endl;
	for (int p = 0; p < Npar; p++) {
		cout << GlobalBestParams.at(p) << "   ";
	}
	cout <<endl<< " --> activations: " << endl;
	for (int f = 0; f < Nfoils; f++) {
		act = 0;
		for (int e = 0; e < Npar; e++) {
			act += GlobalBestParams.at(e) * Resp[f].at(e);
		}
		cout << "foil " << f + 1 << ": " << act << endl;
	}
	delete[] OneBee;

	//write final parameters to extra file, to be read in as histogram:
	FILE *fOutfile2;
	fOutfile2 = fopen("MyOutfile2.txt", "wt");
	for (int p = 0; p < Npar; p++) {
		fprintf(fOutfile2, "%10.3f   ", GlobalBestParams.at(p));
	}
	fclose(fOutfile2);

	// wait for user input so output window doesnt break down:
	cin.get();
}