/*
\file	Generator_DE.cpp
\brief	Evolutionary Aglorithm Generators

\date	2017.6.4
*/

#include <algorithm>
#include <float.h>
#include "Gen.h"


#if defined(WIN32)
#define wxFinite(n) _finite(n)
#elif defined(_LINUX)
#define wxFinite(n) finite(n)
#else
#define wxFinite(n) ((n)==(n))
#endif

namespace az
{
	namespace mea
	{
		namespace gen
		{

			void DE(CPopulationMO& pop, unsigned int idx, CIndividualMO& indnew, double f, double cr)
			{
				unsigned int r1, r2, r3, jrnd;
				do{ r1 = rnd::rand((unsigned int)(0), pop.Size()); } while (r1 == idx);
				do{ r2 = rnd::rand((unsigned int)(0), pop.Size()); } while (r2 == idx || r2 == r1);
				do{ r3 = rnd::rand((unsigned int)(0), pop.Size()); } while (r3 == idx || r3 == r2 || r3 == r1);

				// generate one solution
				jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
				for (unsigned int j = 0; j<pop.P().XSize(); j++)
				{
					if (rnd::rand(0.0, 1.0)< cr || j == jrnd)
						indnew[j] = pop[r1][j] + f*(pop[r2][j] - pop[r3][j]);
					else
						indnew[j] = pop[idx][j];

					if (indnew[j]>pop.P().XUpp(j)) 		indnew[j] = rnd::rand()<0.5 ? 0.5*(pop[idx][j] + pop.P().XUpp(j)) : pop.P().XUpp(j);
					else if (indnew[j]<pop.P().XLow(j))	indnew[j] = rnd::rand()<0.5 ? 0.5*(pop[idx][j] + pop.P().XLow(j)) : pop.P().XLow(j);
				}
			}

			// class XDE
			// constructor
			XDE::XDE()
			{
				mF = 1.0;
				mCR = 1.0;
			}

			// set parameters
			void XDE::Set(double f, double cr)
			{
				mF = f;
				mCR = cr;
			}

			// generate new trial solutions
			CPopulationMO& XDE::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
			{
				unsigned int i, j, r1, r2, r3, jrnd, index;

				//CIndividualMO ind(pop.P());

				popnew.Resize(sizenew);

				//generate new solutions
				for (index = 0; index < sizenew; index++)
				{
					do{ r1 = rnd::rand((unsigned int)(0), pop.Size()); } while (r1 == index);
					do{ r2 = rnd::rand((unsigned int)(0), pop.Size()); } while (r2 == index || r2 == r1);
					do{ r3 = rnd::rand((unsigned int)(0), pop.Size()); } while (r3 == index || r3 == r2 || r3 == r1);

					// generate one solution
					jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
					for (j = 0; j < pop.P().XSize(); j++)
					{
						if (rnd::rand(0.0, 1.0)<mCR || j == jrnd)
							popnew[index][j] = pop[r1][j] + mF*(pop[r2][j] - pop[r3][j]);
						else
							popnew[index][j] = pop[index][j];

						/*if (ind[j]>pop.P().XUpp(j)) 		ind[j] = 0.5*(pop[index][j] + pop.P().XUpp(j));
						else if (ind[j]<pop.P().XLow(j))	ind[j] = 0.5*(pop[index][j] + pop.P().XLow(j));*/

						if (popnew[index][j]>pop.P().XUpp(j))
						{
							double tmp = pop.P().XLow(j) > (2.0*pop.P().XUpp(j) - popnew[index][j]) ? pop.P().XLow(j) : (2.0*pop.P().XUpp(j) - popnew[index][j]);
							popnew[index][j] = rnd::rand() < 0.5 ? tmp : pop.P().XUpp(j);
						}
						else if (popnew[index][j] < pop.P().XLow(j))
						{
							double tmp = pop.P().XUpp(j) < (2.0*pop.P().XLow(j) - popnew[index][j]) ? pop.P().XUpp(j) : (2.0*pop.P().XLow(j) - popnew[index][j]);
							popnew[index][j] = rnd::rand() < 0.5 ? tmp : pop.P().XLow(j);
						}

					}
				}

				return popnew;
			}

			// class XNSDE
			// constructor
			XNSDE::XNSDE()
			{
				mF = 0.8;
				mCR = 0.4;
			}

			// set parameters
			void XNSDE::Set(double f, double cr)
			{
				mF = f;
				mCR = cr;
			}

			// generate new solutions
			CPopulationMO& XNSDE::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
			{
				unsigned int j, r, r1, r2, r3, jrnd, index;

				popnew.Resize(sizenew);

				//generate new solutions
				for (index = 0; index<popnew.Size(); index++)
				{
					//do{r1=rnd::rand((unsigned int)(0), pop.Size());}while(r1==index);
					//do{r2=rnd::rand((unsigned int)(0), pop.Size());}while(r2==index||r2==r1);
					//do{r3=rnd::rand((unsigned int)(0), pop.Size());}while(r3==index||r3==r2||r3==r1);

					do{ r1 = rnd::rand((unsigned int)(0), pop.Size()); r = rnd::rand((unsigned int)(0), pop.Size()); if (pop[r].Rank() < pop[r1].Rank()) r1 = r; } while (r1 == index);
					do{ r2 = rnd::rand((unsigned int)(0), pop.Size()); r = rnd::rand((unsigned int)(0), pop.Size()); if (pop[r].Rank() < pop[r2].Rank()) r2 = r; } while (r2 == index || r2 == r1);
					do{ r3 = rnd::rand((unsigned int)(0), pop.Size()); r = rnd::rand((unsigned int)(0), pop.Size()); if (pop[r].Rank() < pop[r3].Rank()) r3 = r; } while (r3 == index || r3 == r2 || r3 == r1);

					// generate one solution
					jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
					for (j = 0; j<pop.P().XSize(); j++)
					{
						if (rnd::rand(0.0, 1.0)<mCR || j == jrnd)
							popnew[index][j] = pop[r1][j] + mF*(pop[r2][j] - pop[r3][j]);
						else
							popnew[index][j] = pop[index][j];

						// if not a leagle float
						if (wxFinite(popnew[index][j]) == 0) popnew[index][j] = rnd::rand(popnew.P().XLow(j), popnew.P().XUpp(j));

						if (popnew[index][j]>pop.P().XUpp(j)) 		popnew[index][j] = rnd::rand() < 0.5 ? 0.5*(pop[index][j] + pop.P().XUpp(j)) : pop.P().XUpp(j);
						else if (popnew[index][j]<pop.P().XLow(j))	popnew[index][j] = rnd::rand() < 0.5 ? 0.5*(pop[index][j] + pop.P().XLow(j)) : pop.P().XLow(j);

						//if(popnew[index][j]>pop.P().XUpp(j)) 		popnew[index][j] = pop.P().XUpp(j);
						//else if(popnew[index][j]<pop.P().XLow(j))	popnew[index][j] = pop.P().XLow(j);
					}
				}
				return popnew;
			}

			// class XSaDE1
			// constructor
			XSaDE1::XSaDE1()
			{
				mTao1 = 1.0;
				mTao2 = 1.0;
				mTao3 = 1.0;
			}

			// set parameters
			void XSaDE1::Set(double tao1, double tao2, double tao3)
			{
				mTao1 = tao1;
				mTao2 = tao2;
				mTao3 = tao3;
			}

			// generate new solutions
			CPopulationMO& XSaDE1::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop, CPopulationMO& archive)
			{
				CIndividualMO ind(pop.P());

				popnew.Clear();

				std::vector<double> crowdist;
				if (archive.Size() >= 2)
				{
					crowdist.resize(archive.Size());
					for (unsigned int i = 0; i < pop.P().FSize(); i++)
					{
						std::vector<int> index;
						for (unsigned int j = 0; j < archive.Size(); j++) index.push_back(j);

						for (unsigned int j = 0; j < archive.Size()-1; j++)
						for (unsigned int k = j + 1; k < archive.Size(); k++)
						if (archive[index[j]].F(i) > archive[index[k]].F(i))
							std::swap(index[j], index[k]);

						double interval = archive[index[index.size() - 1]].F(i) - archive[index[0]].F(i);
						for (unsigned int j = 1; j<(unsigned int)index.size() - 1; j++) crowdist[index[j]] += (archive[index[j + 1]].F(i) - archive[index[j - 1]].F(i)) / interval;
						crowdist[index[0]] = 1.0E+30;
						crowdist[index[index.size() - 1]] = 1.0E+30;
					}
				}

				//generate new solutions
				for (unsigned int i = 0; i<pop.Size(); i++) pop[i].Rank(0);
				pop.IsSort(false);

				for (unsigned int i = 0; i < sizenew; i++)
				{
					if (pop[i].mOpt == 0)
					{
						// 1: DE/rand/1/bin
						// 2: DE/best/1/bin
						pop[i].mOpt = rnd::rand() <= 0.5 ? 1 : 2;
						//pop[i].mOpt = 1;
					}

					if (pop[i].mAdjust){ AdjustParam(pop[i]); pop[i].mAdjust = false; }

					if (pop[i].mOpt == 1){ DE(pop, i, ind, pop[i].mIf, pop[i].mIcr); }
					else if (pop[i].mOpt == 2) 
					{
						CIndividualMO best(pop.P());

						if (archive.Size() == 0)
						{
							for (unsigned int j = 0; j < pop.P().XSize(); j++) { best[j] = rnd::rand(pop.P().XLow(j), pop.P().XUpp(j)); }
					
						}
						else if (archive.Size() == 1)
						{
							best = archive[0];
						}
						else
						{
							unsigned int b1, b2;
							b1 = rnd::rand((unsigned int)(0), archive.Size());
							do{ b2 = rnd::rand((unsigned int)(0), archive.Size()); } while (b1 == b2);
							best = crowdist[b1] > crowdist[b2] ? archive[b1] : archive[b2];
						}

						DEbest1(pop, i, best, ind, pop[i].mIf, pop[i].mIcr); 
					}

					ind.Evaluate();

					ind.mID = pop[i].mID + sizenew;
					ind.mOpt = pop[i].mOpt;
					ind.mIf = pop[i].mIf;
					ind.mIcr = pop[i].mIcr;

					popnew.Combine(ind);
				}
				return popnew;
			}

			// parameter adjustment
			void XSaDE1::AdjustParam(CIndividualMO& ind)
			{
				//jDE: F_l = 0.1, F_u = 0.9;
				if (rnd::rand() < mTao1)	ind.mIf = 0.1 + rnd::rand()*0.9;
				if (rnd::rand() < mTao2)	ind.mIcr = rnd::rand();
				if (rnd::rand() >= mTao3)	ind.mOpt = rnd::rand() < 0.5 ? 1 : 2;
			}

			// DE/best/1
			void XSaDE1::DEbest1(CPopulationMO& pop, unsigned int idx, CIndividualMO& best, CIndividualMO& indnew, double f, double cr)
			{
				unsigned int r1, r2, jrnd;
				do{ r1 = rnd::rand((unsigned int)(0), pop.Size()); } while (r1 == idx);
				do{ r2 = rnd::rand((unsigned int)(0), pop.Size()); } while (r2 == idx || r2 == r1);

				// generate one solution
				jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
				for (unsigned int j = 0; j<pop.P().XSize(); j++)
				{
					if (rnd::rand(0.0, 1.0)< cr || j == jrnd)
						indnew[j] = best[j] + f*(pop[r1][j] - pop[r2][j]);
					else
						indnew[j] = pop[idx][j];

					if (indnew[j]>pop.P().XUpp(j)) 		indnew[j] = rnd::rand()<0.5 ? 0.5*(pop[idx][j] + pop.P().XUpp(j)) : pop.P().XUpp(j);
					else if (indnew[j]<pop.P().XLow(j))	indnew[j] = rnd::rand()<0.5 ? 0.5*(pop[idx][j] + pop.P().XLow(j)) : pop.P().XLow(j);
				}
			}
		} //namespace gen
	} //namespace mea
} //namespace az
