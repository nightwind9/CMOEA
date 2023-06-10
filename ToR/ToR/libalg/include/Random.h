/*!
\file	Random.h
\brief	random generator

\date	2017.6.4
*/
#define WIN32_LEAN_AND_MEAN

#ifndef	AZ_RANDOM_H
#define	AZ_RANDOM_H

//!\brief	az namespace, the top namespace
namespace az
{

	//!\brief rnd namespace for random generator
	namespace rnd
	{
		//!\brief	initialize the random seed
		//!\param	seeds seeds
		//!\return	void
		void seed(long seeds = 0);

		//!\brief	create a real random in (0.0,1.0)
		//!\return	random number
		double rand();

		//!\brief	create a real random in (low,up) for real number and [low, up) for integer number
		//!\param	low lower range
		//!\param	up upper range
		//!\return	random number
		double rand(double low, double up);

		//!\brief	create a real random in (low,up) for real number and [low, up) for integer number
		//!\param	low lower range
		//!\param	up upper range
		//!\return	random number
		int rand(int low, int up);

		//!\brief	create a real random in (low,up) for real number and [low, up) for integer number
		//!\param	low lower range
		//!\param	up upper range
		//!\return	random number
		unsigned int rand(unsigned int low, unsigned int up);

		//!\brief	create a real Gaussian random with distribution (0.0,1.0)
		//!\return	random number
		double gaussian();

		//!\brief	create a real number with Triangular distribution
		//!\param	min lower limit
		//!\param	max upper limit
		//!\param	mode most likely value 
		//!\return	random number
		double triangular(double min, double max, double mode);

	} //namespace rnd

} //namespace az

#endif //AZ_RANDOM_H
