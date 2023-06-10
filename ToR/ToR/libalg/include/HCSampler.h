/*
\file	HCSampler.h
\brief	samplers in hyper cube

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef	AZ_HCSAMPLER_H
#define	AZ_HCSAMPLER_H

#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief	alg namespace, contains algorithms
	namespace alg
	{
		//!\breif	Latin Hyper Cube design
		//!\param	rand	rand matrix within [low, upp]
		//!\param	low		lower range
		//!\param	upp		upper range
		//!\return	void
		void LHC(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp);

		//!\breif	uniform design
		//!\param	rand	rand matrix within [low, upp]
		//!\param	low		lower range
		//!\param	upp		upper range
		//!\return	void
		void Uniform(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp);

	} //namespace alg

} //namespace az

#endif //AZ_HCSAMPLER_H
