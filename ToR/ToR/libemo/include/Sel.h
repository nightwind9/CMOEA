/*
\file	Sel.h
\brief	MOEA selection strategies

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef AZ_SELECTION_H
#define AZ_SELECTION_H

#include <vector>
#include <list>
#include "IndividualMO.h"
#include "PopulationMO.h"
//#include "Initialization.h"
//#include "alg/amoeba.h"

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{
		//!\brief select strategies
		namespace sel
		{

			//!\brief an empty selection (noting to do)
			class SEmpty
			{
			public:
				//!\brief	select from current population and offspring population
				//!\param	pop combined population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Select(CPopulationMO& pop, unsigned int size) { return pop; }
			};//class SEmpty

			//!\brief original crowded selection strategy in NSGA-II
			class SCrowd
			{
			public:
				//!\brief	sort population, the first 'size' are best ones
				//!\param	pop population
				//!\param	size size of best ones
				//!\return	population
				CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

				//!\brief	select from current population and offspring population
				//!\param	pop combined population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
			};//class SCrowd

			//////////////////////////////////////////////////////////////////////////////////////////////////////
			//!\brief MODIFIED crowded selection strategy in NSGA-II
			class SCrowd2
			{
			public:
				//!\brief	sort population, the first 'size' are best ones
				//!\param	pop population
				//!\param	size size of best ones
				//!\return	population
				CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

				//!\brief	select from current population and offspring population
				//!\param	pop combined population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
			};//class SCrowd2


			//!\ Double crowding distance ranking (for ToR) 
			//!\ R1: without considering constraints.
			//!\ R2: CDP
			class SDouCrowd
			{
			public:
				//!\brief	sort population, the first 'size' are best ones
				//!\param	pop population
				//!\param	size size of best ones
				//!\return	population
				CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

				//!\brief	select from current population and offspring population
				//!\param	pop combined population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Select(CPopulationMO& pop, unsigned int size);

			private:
				//!\brief	assign rank value and sort population
				//!\param	pop population
				//!\param	conrank indicate that whether the conducting rank considers constraints
				//!\return	vector 
				std::vector<unsigned int> DouRank(CPopulationMO& pop, bool conrank);

				//!\brief	calculate crowding distance values of solutions in the same rank
				//!\param	pop		current population
				//!\param	rk		current rank obtained by non-dominated sort
				//!\reutrn	void
				void CrowdDistSort(CPopulationMO& pop, std::vector<unsigned int>& rk);

			}; // class SDouCrowd

		}//namespace sel

	} //namespace mea

} //namespace az


#endif //AZ_SELECTION_H
