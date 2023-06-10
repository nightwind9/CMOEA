/*
\file	IndividualMO.h
\brief	individual class for MOEA

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef	AZ_INDIVIDUALMO_H
#define	AZ_INDIVIDUALMO_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Parameter.h"
#include "Random.h"

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{

		//!\brief individual class for MOEA
		class CIndividualMO
		{
		public:
			unsigned int mRank,			//!< rank value
				mID,					//!< ID
				mOpt;					//!< ID of Optimizer, used for hybrid method
			double	mConstraint;		//!< constraint = sum(|mvEq|) + sum(mvIneq > 0)

			std::vector<double>	mvX,	//!< variables
				mvF,					//!< objective values
				mvEq,					//!< equality values
				mvIneq;					//!< inequality values
			CParameter* pPar;			//!< the pointer to the parameter

			double  mIf, mIcr;			//!< F and CR for self-adaptive DE
			bool	mAdjust;			//!< indicate whether the DE parameters need to adjust

		public:

			//!\brief	Constractor
			//!\param	par parameters
			//!\return	void
			CIndividualMO(CParameter& par);

			//!\brief	Constractor
			//!\param	ind another individual
			//!\return	void
			CIndividualMO(CIndividualMO& ind);

			//!\brief	get the parameter object reference
			//!\return	the parameter object reference
			inline CParameter& P() { return *pPar; }

			//!\brief	check to see if it's feasible
			//!\return	whether the indiviudal is feasible
			inline bool IsFeasible() { return mConstraint <= P().TolC(); }

			//!\brief	check to see if it's epsilon-feasible
			//!\return	whether the indiviudal is epsilon-feasible
			inline bool IsEFeasible() { return mConstraint <= P().EPSN_CHT() + P().TolC(); }
			//inline bool IsEFeasible() { return mConstraint <= P().TolC(); }

			//!\brief	get the ith objective
			//!\param	i objective index
			//!\return	the ith objective reference
			inline double& F(unsigned int i) { return mvF[i]; }

			//!\brief	get the ith variable
			//!\param	i variable index
			//!\return	the ith variable reference
			inline double& X(unsigned int i) { return mvX[i]; }


			//!\brief	get the ith equality constraint
			//!\param	i objective index
			//!\return	the ith objective reference
			inline double& E(unsigned int i) { return mvEq[i]; }

			//!\brief	get the ith inequality constraint
			//!\param	i variable index
			//!\return	the ith variable reference
			inline double& I(unsigned int i) { return mvIneq[i]; }

			//!\brief	get the ith variable
			//!\param	i variable index
			//!\return	the ith variable reference
			inline double& operator[](unsigned int i) { return mvX[i]; }

			//!\brief	get the constraint reference
			//!\return	the constraint reference
			inline double& C() { return mConstraint; }

			//!\brief	get the objective vector reference
			//!\return	the objective vector reference
			inline std::vector<double>&	F() { return mvF; }

			//!\brief	get the variable vector reference
			//!\return	the variable vector reference
			inline std::vector<double>&	X() { return mvX; }

			//!\brief	get the equality constraint vector reference
			//!\return	the objective vector reference
			inline std::vector<double>&	E() { return mvEq; }

			//!\brief	get the inequality constraint vector reference
			//!\return	the variable vector reference
			inline std::vector<double>&	I() { return mvIneq; }

			//!\brief	get the variable vector reference
			//!\return	the variable vector reference
			inline std::vector<double>&	operator()() { return mvX; }

			//!\brief	set the rank value
			//!\param	r new rank value
			//!\return	rank value
			inline unsigned int Rank(unsigned int r) { mRank = r; return r; }

			//!\brief	get the rank value
			//!\return	rank value
			inline unsigned int Rank() { return mRank; }

			//!\brief	set ID
			//!\param	id new ID
			//!\return	new ID
			inline unsigned int ID(unsigned int id) { mID = id; return id; }

			//!\brief	get ID
			//!\return	ID
			inline unsigned int ID() { return mID; }

			//!\brief	set Optizer
			//!\param	opt new Optimizer
			//!\return	new Optimizer
			inline unsigned int OPT(unsigned int opt) { mOpt = opt; return opt; }

			//!\brief	get Optimizer
			//!\return	Optimizer
			inline unsigned int OPT() { return mOpt; }

			//!\brief	dominance check
			//!\param	ind another individual
			//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
			int Dominate(CIndividualMO& ind);

			//!\brief	constratint-dominance check
			//!\param	ind another individual
			//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
			int CDominate(CIndividualMO& ind);

			//!\brief	e-constratint-dominance check
			//!\param	ind another individual
			//!\param	e value of epsilon
			//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
			int CEDominate(CIndividualMO& ind);

			//!\brief	evaluate the individual
			//!\return	void
			void Evaluate();

			//!\brief	make indiviudal to be a feasible one
			//!\brief	no
			void Check();

			//!\brief	set variable to be another individual
			//!\param	ind another individual
			//!\return	reference to the individual
			CIndividualMO& operator =(CIndividualMO& ind);

			//!\brief	check to see if two individuals are equal
			//!\param	ind another individual
			//!\return	whether they are equal
			bool operator ==(CIndividualMO& ind);

			//!\brief	check to see which is better
			//!\param	ind another individual
			//!\return	true if this individual is better than ind
			bool operator<(CIndividualMO& ind);

			//!\brief	write to a stream
			//!\param	os output stream
			//!\param	ind individual
			//!\return	output stream
			friend std::ostream& operator<<(std::ostream& os, CIndividualMO& ind);

			//!\brief	read from a stream
			//!\param	is input stream
			//!\param	ind individual
			//!\return	input stream
			friend std::istream& operator>>(std::istream& is, CIndividualMO& ind);
		};//class CIndividualMO

	} //namespace mea

} //namespace az

#endif //AZ_INDIVIDUALMO_H
