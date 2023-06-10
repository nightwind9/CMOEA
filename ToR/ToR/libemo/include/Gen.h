/*
\file	Gen.h
\brief	Evolutionary Aglorithm Generators

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef AZ_GENERATOR_EA_H
#define AZ_GENERATOR_EA_H

#include "IndividualMO.h"
#include "PopulationMO.h"
#include "Initialization.h"

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{

		//!\brief	gen namespace, offspring generate strategies
		namespace gen
		{

			//!\brief	Polynomial mutation operator (PM)
			//\param	ind individual to be mutated
			//\return	indivdual
			CIndividualMO& PM(CIndividualMO& ind);

			//!\brief	Simulated Binary crossover (SBX)
			//\param	son1 offspring
			//\param	son2 offspring
			//\param	parent1 parent
			//\param	parent2 parent
			//\return	void
			void SBX(CIndividualMO& son1, CIndividualMO& son2, CIndividualMO& parent1, CIndividualMO& parent2);

			//!\brief	Differential Evolution (DE)
			//!\param	pop	current population
			//!\param	idx	current individual to be evoluted
			//!\param	indnew	generated individual
			//!\param	f	scalar factor
			//!\param	cr crossover probability
			//!\return	void
			void DE(CPopulationMO& pop, unsigned int idx, CIndividualMO& indnew, double f, double cr);

			/////////////////////////////////////////////////////////////////////////////////////
			//!\brief	DE operator
			class XDE
			{
			protected:
				double	mF,				//!< setp length
					mCR;			//!< crossover probability
			public:
				//!\brief	constructor
				//!\brief	no
				XDE();

				//!\brief	set parameters
				//!\param	f		 step length
				//!\param	cr		 crossover probability
				//!\return	void
				void Set(double f, double cr);

				//!\brief	Generator
				//!\param	sizenew size of offpsring population
				//!\param	popnew	offspring population
				//!\param	pop		parent population
				//!\return	offspring population
				CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
			};//class DE

			/////////////////////////////////////////////////////////////////////////////////////////////
			//!\brief	SBX
			class XSBX
			{
			public:
				//!\brief	Generator
				//!\param	sizenew size of offpsring population
				//!\param	popnew	offspring population
				//!\param	pop		parent population
				//!\return	offspring population
				CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
			};//class SBX

			/////////////////////////////////////////////////////////////////////////////////////
			//!\brief	NSDE operator
			class XNSDE
			{
			protected:
				double	mF,				//!< setp length
						mCR;			//!< crossover probability
			public:
				//!\brief	constructor
				//!\brief	no
				XNSDE();

				//!\brief	set parameters
				//!\param	f		 step length
				//!\param	cr		 crossover probability
				//!\return	void
				void Set(double f, double cr);

				//!\brief	Generator
				//!\param	sizenew size of offpsring population
				//!\param	popnew	offspring population
				//!\param	pop		parent population
				//!\return	offspring population
				CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
			};//class XNSDE

			/////////////////////////////////////////////////////////////////////////////////////
			//!\brief	Self-adaptive DE operator
			class XSaDE1
			{
			protected:
				double mTao1, mTao2, mTao3;

			public:
				//!\brief	constructor
				//!\brief	no
				XSaDE1();

				//!\brief	set parameters
				//!\param	tao1	probability to adjust F 
				//!\param	tao2	probability to adjust CR
				//!\param	tao3	probability to aDjust trial vector generation strategy
				//!\return	void
				void Set(double tao1, double tao2, double tao3);

				//!\brief	Generator
				//!\param	sizenew size of offpsring population
				//!\param	popnew	offspring population
				//!\param	pop		parent population
				//!\return	offspring population
				CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop, CPopulationMO& archive);

			protected:
				//!\brief	DE/best/1
				//!\param	pop	current population
				//!\param	idx	current individual to be evoluted
				//!\param	best the best individual found
				//!\param	indnew	generated individual
				//!\param	f	scalar factor
				//!\param	cr crossover probability
				//!\return	void
				void DEbest1(CPopulationMO& pop, unsigned int idx, CIndividualMO& best, CIndividualMO& indnew, double f, double cr);

				//!\brief	self-adaptive adjustment of parameters
				//!\param	ind	current individual
				//!\return	void
				void AdjustParam(CIndividualMO& ind);

			};//class DE

		} //namespace gen

	} //namespace mea

} //namespace az

#endif //AZ_GENERATOR_EA_H
