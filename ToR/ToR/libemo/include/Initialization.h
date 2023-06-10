/*
\file	Initialization.h
\brief	initialization strategies

\author Zhongwei Ma

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef AZ_INITIALIZATION_H
#define AZ_INITIALIZATION_H

#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{
		//!\brief initialization strategies
		namespace ini
		{
			//!\brief hybrid initialization strategy
			class Hybrid
			{
			protected:
				unsigned int		mEvas;		//!< evaluations
			public:
				Hybrid();

				~Hybrid();

				//!\brief	calculate the scalar obj
				//!\brief	obj value
				double obj(double* xy, int Dimension);

				//!\brief	get evaluations
				//!\return	evaluations
				inline unsigned int EvaTimes() { return mEvas; }

				//!\brief	initialize a population
				//!\param	pop population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);
			};//class Hybrid

			//!\brief random initialization strategy
			class Uniform
			{
			protected:
				unsigned int mEvas;	//!< evaluations
			public:
				//!\brief	get evaluations
				//!\return	evaluations
				inline unsigned int EvaTimes() { return mEvas; }

				//!\brief	initialize a population
				//!\param	pop population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);
			};//class Random

			//!\brief initialization with experiment design
			class LHC
			{
			protected:
				unsigned int mEvas;	//!< evaluations
			public:
				//!\brief	get evaluations
				//!\return	evaluations
				inline unsigned int EvaTimes() { return mEvas; }

				//!\brief	initialize a population
				//!\param	pop population
				//!\param	size population size
				//!\return	population
				CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);
			};//class RandomSet

		}//namespace ini

	} //namespace mea

} //namespace az


#endif //AZ_INITIALIZATION_H
