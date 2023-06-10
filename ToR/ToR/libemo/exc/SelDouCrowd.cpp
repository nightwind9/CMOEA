/*
\file	Selection_Crowd.cpp
\brief	MOEA selection strategies

\date	2017.9.18
*/

#include <vector>
#include <list>
#include "Sel.h"
//#include "emo/Nondominate.h"
//#include "emo/Pareto.h"

namespace az
{
	namespace mea
	{
		namespace sel
		{
			// Scrowd1:
			// select from current population and offspring population
			CPopulationMO& SDouCrowd::Select(CPopulationMO& pop, unsigned int size)
			{
				if (pop.Size() > size) SelectSort(pop, size).Erase(size); return pop;
			}

			// sort population, the first 'size' are best ones
			CPopulationMO& SDouCrowd::SelectSort(CPopulationMO& pop, unsigned int size)
			{
				if (pop.Size() <= size) return pop;

				unsigned int start, end;
				std::vector<unsigned int> R1, R2, R3;

				// Step 1: rank sort
				// CDP
				bool crank = true;
				R1 = DouRank(pop, crank);
				// without considering constraints
				crank = false;
				R2 = DouRank(pop, crank);
				
				// Step 2: calculate the rank R3
				double fea_num = 0.0;
				for (unsigned int i = 0; i < pop.Size(); i++) if (pop[i].IsFeasible()) fea_num += 1.0;
				fea_num = fea_num / (double)(pop.Size());

				for (unsigned int i = 0; i < pop.Size(); i++) R3.push_back(0);
				for (unsigned int i = 0; i < pop.Size(); i++)
				{

					R3[R1[i]] += (0.5 + 0.5*fea_num)*i;
					R3[R2[i]] += (0.5 - 0.5*fea_num)*i;
				}

				// Step 3: sort the population
				for (unsigned int i = 0; i < pop.Size()-1; i++)
				for (unsigned int j = pop.Size() - 1; j > i; j--)
				{
					if (R3[j] < R3[j - 1])
					{
						std::swap(R3[j - 1], R3[j]);
						pop.Swap(j - 1, j);
					}
				}

				pop.IsSort(true);
				return pop;
			}

			// assign rank value and sort population
			std::vector<unsigned int> SDouCrowd::DouRank(CPopulationMO& pop, bool conrank)
			{
				std::vector<unsigned int> R;
				for (unsigned int i = 0; i < pop.Size(); i++) { R.push_back(i); }

				// has been sorted before
				if (pop.IsSort()) return R;

				// only has one individual
				if (pop.Size() < 2) { pop.IsSort(true); return R; }

				int s, t, better;
				unsigned int minRank, curRank = 0, noAssign, size;

				s = 0; t = pop.Size() - 1;
				size = pop.Size();

				/*if (conrank)
				{
					while (s<t)
					{
						while (s < int(pop.Size()) && pop.In(R[s]).IsFeasible()) s++;
						while (t >= 0 && !pop.In(R[t]).IsFeasible()) t--;
						if (s<t)	{ std::swap(R[s], R[t]); s++; t--; }
					}
					size = s;
				}*/

				std::vector<bool>						vAssign(size);
				std::vector<unsigned int>				v2Assign;
				std::vector<unsigned int>				vBeDom(size);
				std::vector<std::vector<unsigned int>>	vDom(size);
				std::vector<unsigned int>::iterator		it1, it2;

				//Step 2: Initialize
				for (unsigned int i = 0; i < size; i++) vAssign[R[i]] = false;

				// Step 3: MOGAFonseca flow
				for (unsigned int i = 0; i < size; i++)
				{
					for (unsigned int j = i + 1; j < size; j++)
					{
						better = conrank ? pop.In(R[i]).CDominate(pop.In(R[j])) : pop.In(R[i]).Dominate(pop.In(R[j]));
						//R[i] is dominated by R[j]
						if (better<0) { vBeDom[R[i]]++; vDom[R[j]].push_back(R[i]); }
						//R[j] is dominated by R[i]
						else if (better>0) { vBeDom[R[j]]++; vDom[R[i]].push_back(R[j]); }
					}	//end for
				}	// end for

				// Step 4: Assign rank to feasible solutions
				noAssign = size;
				while (noAssign > 0)
				{
					curRank++;

					minRank = size;
					//find the cluster to assign a rank value
					for (unsigned int i = 0; i<size; i++) if (!vAssign[R[i]])
					{
						if (vBeDom[R[i]]<minRank)
						{
							minRank = vBeDom[R[i]];
							v2Assign.clear(); v2Assign.push_back(R[i]);
						}
						else if (vBeDom[R[i]] == minRank) v2Assign.push_back(R[i]);
					}

					//CHECK( v2Assign.size()>0, "CRankSort::Sort()" );

					//assign rank
					it1 = v2Assign.begin();
					while (it1 != v2Assign.end())
					{
						vAssign[*it1] = true;
						pop.In(*it1).Rank(curRank);
						it2 = vDom[*it1].begin();
						while (it2 != vDom[*it1].end()) vBeDom[*it2++]--;
						it1++; noAssign--;
					}//for
				}//end while

				//Step 4: Sort the population by rank
				for (s = 0; s < size - 1; s++)
				for (t = s + 1; t < size; t++)
				if (pop.In(R[t]).Rank() < pop.In(R[s]).Rank())
				{
					std::swap(R[s], R[t]);
				}

				/*if(conrank) for (unsigned int i = 0; i < pop.Size(); i++)
				{
					if (pop.In(i).IsFeasible()) curRank = pop.In(i).Rank();
					else pop.In(i).Rank(curRank + 1);
				}*/
				
				CrowdDistSort(pop, R);

				return R;
			}// SDouCrowd()

			// calculate crowding distance values of solutions in the same rank
			void SDouCrowd::CrowdDistSort(CPopulationMO& pop, std::vector<unsigned int>& rk)
			{
				unsigned int start = 0, end = 0;
				while (end < pop.Size())
				{
					start = end;
					while (end<pop.Size() && pop[rk[end]].Rank() == pop[rk[start]].Rank()) end++;

					if (end > start + 2 && start < pop.Size() - 2)
					{
						unsigned int i, j, k;
						double interval;

						std::vector<double>			share(end - start);
						std::vector<unsigned int>		index(end - start);

						for (i = 0; i<(unsigned int)share.size(); i++) share[i] = 0.0;

						//calculate the share values
						for (i = 0; i<pop.P().FSize(); i++)
						{
							for (j = start; j<end; j++) index[j - start] = j;

							for (j = 0; j<(unsigned int)index.size() - 1; j++)
							for (k = j + 1; k<(unsigned int)index.size(); k++)
							if (pop[rk[index[j]]].F(i) > pop[rk[index[k]]].F(i))
								std::swap(index[j], index[k]);

							interval = pop[rk[index[index.size() - 1]]].F(i) - pop[rk[index[0]]].F(i);
							for (j = 1; j<(unsigned int)index.size() - 1; j++) share[index[j] - start] += (pop[rk[index[j + 1]]].F(i) - pop[rk[index[j - 1]]].F(i)) / interval;
							share[index[0] - start] = MAXDOUBLE;
							share[index[index.size() - 1] - start] = MAXDOUBLE;
						}

						//sort the sub-population according to the share value
						for (i = start; i<end; i++)
						for (j = i + 1; j<end; j++)
						if (share[i - start] < share[j - start])
						{
							std::swap(rk[i], rk[j]);
							std::swap(share[i - start], share[j - start]);
						}

						share.clear();
						index.clear();
					}// end if
				}// end while
			}// CrowdDistSort1()

		}//namespace sel
	} //namespace mea
} //namespace az
