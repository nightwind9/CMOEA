/*
\file	Metric.h

\date	2017.10.18
*/
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <cassert>
#include <iomanip>
using namespace std;

//!\brief	az namespace, the top namespace
namespace az
{
	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{
		//!\brief namespace of constrained evolutionary algoirhtm 
		namespace cmea
		{
			// class EVAL
			class EVAL
			{
			protected:
				vector<vector<double>>  pf;			//!< vector to store the true PF
				vector<vector<double>>  result;		//!< vector to store the population to be evaluated

				unsigned int eMaxRun,				//!< number of runs
					eMaxStep,						//!< maximal generation
					ePopSize,						//!< population size
					eObjNum;						//!< number of objectives
				double igd_value,					//!< value of IGD metric 
					gd_value,						//!< value of GD metric
					sp_value,						//!< value of SP metric
					ms_value,						//!< value of MS metric
					hv_value;						//!< value of HV metric
				vector<vector<double>> eTrackIGD,		//!< evolution of IGD metric
					eFeaRatio;						//!< evolution of feasible ratio

				string eProbName;					//!< problem name

			public:
				//!\brief	constractor
				//!\param	maxrun		the number of parallel implements
				//!\param	maxstep		maximal generation
				//!\param	popsize		population size
				//!\param	objnum		number of objectives
				//!\param	probname	problem name
				EVAL(unsigned int	maxrun,
					unsigned int	maxstep,
					unsigned int	popsize,
					unsigned int	objnum,
					string			probname);

				//!\brief	calculate the IGD value of each generation
				//!\param	CurRun		index of current run
				//!\param	CurStep		index of current step
				//!\param	pop			population to be evaluated
				//!\param	feasible	the feasible ration of current popualtion
				//!\return  none
				void EvolutionIGD(unsigned int CurRun, unsigned int CurStep, vector<vector<double>> pop, double fearatio);

				//!\brief	evaluate the final popualtion
				//!\param	probname	problem name
				//!\return  none
				void FinalResults(vector<vector<double>> pop);

				//!\brief	output the total result
				//!\param	probname	problem name
				//!\return  none
				void OutputTotal();

			protected:
				//!\brief	load true PF from txt files
				//!\param	probname	problem name
				//!\param	pf			vector to store data
				//!\param	nobj		number of objectives
				//!\return	none
				void loadpfront(string probname, vector< vector<double>> &pf, int nobj);

				//!\brief	calculate the IGD value of a given population
				//!\param	pop		current population
				//!\return	none
				void Inverted_GD(vector<vector<double>> pop);

				//!\brief	calculate the GD value of a given population
				//!\param	pop		current population
				//!\return	none
				void Gen_Distance(vector<vector<double>> pop);

				//!\brief	calculate the SP value of a given population
				//!\param	pop		current population
				//!\return	none
				void Spacing(vector<vector<double>> pop);

				//!\brief	calculate the SP value of a given population
				//!\param	pop		current population
				//!\return	none
				void Max_spread(vector<vector<double>> pop);

				void Hypervolume(vector<vector<double>> pop);

				void Hypervolume_slice(vector<pair<double, vector<vector<double>>>> &stemp, vector<vector<double>> p1, int k, vector<double> refpoint);

				void Hypervolume_Insert(vector<double> p, int k, vector<vector<double>> &p1);

				void Hypervolume_Add(pair<double, vector<vector<double>>> cell_, vector<pair<double, vector<vector<double>>>> &S);

				//!\brief	calculate the distance of two vectors
				//!\param	vec1	the reference of the first vector
				//!\param	vec2	the reference of the second vector
				//!\return	the distance
				double dist_vector(vector <double> &vec1, vector <double> &vec2);

				//!\brief	calculate the standard deviation of a given set
				//!\param	vec		the given set 
				//!\return	standard deviation
				double STDev(vector<double> &vec);

				//!\brief	calculate the average value of a given set
				//!\param	vec		the given set 
				//!\return	the average value
				double Avg(vector<double> &vec);


			};

		} //namespace cmea
	} //namespace mea
} //namespace az


