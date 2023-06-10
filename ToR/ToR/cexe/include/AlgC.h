/*
\file	AlgD.h
\brief	Framwork for MOEA

\date	2017.6.5
*/

#pragma once

#include <ctime>
#include <vector>
#include <string>
//#include "alg/Matrix.h"
//#include "alg/GTM.h"
#include "Parameter.h"
#include "PopulationMO.h"
#pragma comment(lib,"libemo.lib")
#pragma comment(lib,"libalg.lib")

#include"engine.h"
#include"memory.h"
#pragma comment(lib,"libmat.lib")
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libeng.lib")

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{

		//!\brief namespace of constrained evolutionary algoirhtm 
		namespace cmea
		{

			//!\brief	Contrained MOEA
			class CMOO
			{
			protected:
				// parameters for algorithm
				unsigned int
					mPopSize,				//!< population size
					mMaxStep,				//!< maximum generation
					mStep,					//!< current step
					mEvas;					//!< the calculation number of the objectives

				double mCP;					//!< the parameter to control the reducing speed of epsilon level

				CPopulationMO	mPop;		//!< current population
				CPopulationMO	ex_pop;		//!< external population
				CParameter*	pPar;			//!< pointer to the parameter object
				std::string	mOptimizer;		//!< optimizer inside the time windows

			public:
				//!\brief	constractor
				//!\param	strategy	prediction strategy
				//!\param	popsize		population size
				//!\param	stepmax		maximum steps
				//!\param	taot		the length of each state in case of generations
				//!\param	nt			servety of changes
				//!\param	torder		order of time series
				//!\parame	t0			initial state of time T
				//!\param	alpha		variance control parameter
				//!\param	par			parameter object
				//!\return	void
				CMOO(
					std::string&	optimizer,
					unsigned int	popsize,
					unsigned int	stepmax,
					CParameter&		par);

				//!\brief	get the pointer to the parameter object
				//!\return	pointer to the parameter object
				inline CParameter& P() { return *pPar; }

				//!\brief	get the population
				//1\return	reference to population
				inline CPopulationMO& Population() { return mPop; }

				//!\brief	get the objective evaluation times
				//!\return	objective evaluation times
				inline unsigned int& EvaTimes() { return mEvas; }

				//!\brief	check to see whether the terminate condition is satified
				//!\return	true if terminate
				inline bool IsTerminate() { return mStep >= mMaxStep; }

				//!\brief	get the current step
				//!\return	current step
				inline unsigned int CurStep() { return mStep; }

				//!\brief	save population
				//!\return	void
				inline void Write(std::string name) { Population().Write(name); }

				//!\brief	initialization
				//!\return	void
				void Init();

				//!\brief	reset parameters
				//!\return	void
				void Reset();

				//!\brief	one step 
				//!\return	current step
				unsigned int Step();

				//!\brief	update external population
				//!\return	void
				void update_external_pop(CPopulationMO& pop, CPopulationMO& ex_pop);

			protected:
				//!\brief	make all new solutions in feasible region
				//!\param	pop offspring population
				//!\return	void
				void Check(CPopulationMO& pop);

				//!\brief	make all new solutions are different from old solutions
				//!\param	popnew offspring population
				//!\param	pop old population
				//!\return	void
				void Check(CPopulationMO& popnew, CPopulationMO& pop);

				//!\brief	generate offsprings
				//!\param	popnew offspring population
				//!\param	size offspring population size
				//!\return	offspring population
				CPopulationMO& Generate(CPopulationMO& popnew, unsigned int size);

				//!\brief	environmental selection
				//!\param	pop current population
				//!\param	size current population size
				//!\return	offspring population
				CPopulationMO& Select(CPopulationMO& pop, unsigned int size);


			}; //class CMOO

			//!\brief	Plot function
			void ePlot(Engine*& ep, std::string& problem, unsigned int obj);

			//!\brief	Print Copyright Information
			void PrintCopyright();

			//!\brief	Print parameter settings
			void ParaSetting(std::string method, unsigned int runs, unsigned int maxgen, unsigned int popsize, unsigned int varsize);

			//!\brief	Show the Rate of Progress
			void ShowProgress(unsigned int run, unsigned int gen, unsigned int ic, unsigned int maxgen);

		} //namespace cmea

	} //namespace mea

} //namespace az
