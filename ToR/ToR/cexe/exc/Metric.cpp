#include "Metric.h"
#include <algorithm>
#include <vector>
#include <random>
#include <iomanip>
#include <iostream>

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
			// constructor
			EVAL::EVAL(unsigned int maxrun,
				unsigned int maxstep,
				unsigned int popsize,
				unsigned int objnum,
				string probname)
			{
				eMaxRun = maxrun;
				eMaxStep = maxstep;
				ePopSize = popsize;
				eObjNum = objnum;
				eProbName = probname;

				// load the true PF of current problem
				loadpfront(eProbName, pf, eObjNum);

				// for (unsigned int i = 0; i < eMaxRun; i++) { vector<double> EachRun(eMaxStep, -1.0); eTrackIGD.push_back(EachRun); }

				// for (unsigned int i = 0; i < eMaxRun; i++) { vector<double> EachRun(maxstep, 0.0); eFeaRatio.push_back(EachRun); }
			}

			// calculate the IGD value of each generation
			void EVAL::EvolutionIGD(unsigned int CurRun, unsigned int CurStep, vector<vector<double>> pop, double fearatio)
			{
				if (eFeaRatio.size() > 0 && CurStep < eMaxStep && CurRun < eMaxRun)	eFeaRatio[CurRun][CurStep] = fearatio;
				Inverted_GD(pop);
				//Gen_Distance(pop);
				if (eTrackIGD.size() > 0 && CurStep < eMaxStep && CurRun < eMaxRun)	eTrackIGD[CurRun][CurStep] = igd_value;
			}

			// evaluate the final population
			void EVAL::FinalResults(vector<vector<double>> pop)
			{
				vector<double> r;

				// Spacing(pop);	r.push_back(sp_value);
				Max_spread(pop);	r.push_back(ms_value);
				Gen_Distance(pop);	r.push_back(gd_value);
				Inverted_GD(pop);	r.push_back(igd_value);

				result.push_back(r);

				r.clear();
			}

			// output the final result
			void EVAL::OutputTotal()
			{
				std::stringstream str;
				str << "Results/Final_" << eProbName.c_str() << ".dat";

				ofstream fout;

				fout.open(str.str().c_str(), ios::app);
				fout << setprecision(8) << setiosflags(ios::fixed);
				if (fout.is_open())
				{
					/*fout << "Track_IGD:"<<"\t"<<"FeaRatio:\t" << endl;
					for (unsigned int i = 0; i < eMaxStep; i++){ fout << eTrackIGD[i] / (double)(eMaxRun) << "\t" << eFeaRatio[i] / (double)(result.size()) << endl; }
					fout << endl;*/

					//fout << "Track_IGD:"<< endl;
					//for (unsigned int i = 0; i < eMaxStep; i++)
					//for (unsigned int j = 0; j < eMaxRun; j++)
					//{
					//	fout << eTrackIGD[j][i];
					//	if (j < eMaxRun - 1)
					//	{
					//		fout << "\t";
					//	}
					//	else
					//	{
					//		fout << endl;
					//	}
					//}
					//fout << endl;

					/*fout << "Statistacis:" << endl;*/
					if (result.size() == 0) std::cout << "  NO feasible soltuions have been found over all runs!" << std::endl;

					if (result.size() != 0)
					{
						std::string m[] = { "MS: ", "GD: ", "IGD: " };
						for (unsigned int i = 0; i < result[0].size(); i++)
						{
							std::vector<double> val;
							for (unsigned int j = 0; j < result.size(); j++) val.push_back(result[j][i]);
							fout << setprecision(3) << scientific << Avg(val) << "(¡À" << STDev(val) << ")" << "	";

							std::cout<< "  " << m[i] << setprecision(2) << scientific << Avg(val) << "(¡À" << STDev(val) << ")" << " ";
							
							//fout << Avg(val) << "(¡À" << STDev(val) << ")" << "	";
						}
					}
					fout << eMaxRun - result.size() << endl;
					std::cout << endl << endl;

					//fout << "SP:	" << "	" << "MS:	" << " " << "GD:	" << "	" << "IGD:	" << endl;
					fout << "MS:\t" << "\t" << "GD:\t" << "\t" << "IGD:\t" << endl;
					if (result.size() != 0) for (unsigned int i = 0; i < result.size(); i++)
					{
						for (unsigned int j = 0; j < result[i].size(); j++)
						{
							fout << result[i][j];
							if (j != result[i].size() - 1)	fout << "	";
						}
						fout << endl;
					}
				}
				else
				{
					std::cout << "failed to open " << str.str().c_str() << endl;
				}
				fout.close();
			}

			//load true PF from txt files
			void EVAL::loadpfront(string probname, vector< vector<double>> &pf, int nobj)
			{
				std::stringstream str0;
				str0 << "TruePF/" << probname.c_str() << ".dat";

				std::fstream fin;
				int line = 0;
				char str[100] = " ";
				fin.open(str0.str().c_str(), std::ios::in);
				if (fin.is_open())
				{
					const char* pch2;
					char  a[20], b[20], c[20], d[20];
					std::string str;

					while (!fin.eof())
					{
						vector<double> data;
						std::getline(fin, str, '\n');
						pch2 = str.c_str();
						sscanf(pch2, "%s %s %s %s", a, b, c, d);
						data.push_back(atof(a));
						data.push_back(atof(b));
						if (nobj == 3)
							data.push_back(atof(c));
						//if(nobj==4) 
						//{
						//	data.y_obj[3] = atof(d);
						//}
						line++;
						pf.push_back(data);
					}
				} //end if
				else
					std::cout << "failed to open " << str0.str().c_str() << endl;
				fin.close();
			}

			// calculate the IGD value of a given population
			void EVAL::Inverted_GD(vector<vector<double>> pop)
			{
				igd_value = 0.0;
				for (int i = 0; i<pf.size(); i++)
				{
					double min_d = 1.0E+30;
					for (int j = 0; j < pop.size(); j++)
					{
						double d = dist_vector(pf[i], pop[j]);
						if (d<min_d)  min_d = d;
					}
					igd_value += min_d;
				}
				igd_value /= pf.size();
			}

			// calculate the GD value of a given population
			void EVAL::Gen_Distance(vector<vector<double>> pop)
			{
				gd_value = 0.0;
				for (int i = 0; i<pop.size(); i++)
				{
					double min_d = 1.0E+30;
					for (int j = 0; j<pf.size(); j++)
					{
						double d = dist_vector(pop[i], pf[j]);
						if (d<min_d)  min_d = d;
					}
					gd_value += min_d;
				}
				gd_value /= pop.size();
			}

			// calculate the SP value of a given population
			void EVAL::Spacing(vector<vector<double>> pop)
			{
				int psize = pop.size();
				vector<double> dist(psize, 1.0E+30);
				for (int i = 0; i < psize; i++)
				{
					for (int j = 0; j < psize; j++)
					{
						if (j != i)
						{
							double tmp = dist_vector(pop[i], pf[j]);
							if (tmp < dist[i]) dist[i] = tmp;
						}
					}
				}
				sp_value = STDev(dist);
			}

			// calculate the SP value of a given population
			void EVAL::Max_spread(vector<vector<double>> pop)
			{
				int pfsize = pf.size();
				int psize = pop.size();
				vector<double> pf_max(eObjNum, -1.0e+30), pf_min(eObjNum, 1.0e+30);
				vector<double> p_max(eObjNum, -1.0e+30), p_min(eObjNum, 1.0e+30);
				for (int n = 0; n < eObjNum; n++)
				{
					for (int i = 0; i < pf.size(); i++)
					{
						if (pf_max[n] < pf[i][n]) pf_max[n] = pf[i][n];
						if (pf_min[n] > pf[i][n]) pf_min[n] = pf[i][n];
					}
					for (int i = 0; i < pop.size(); i++)
					{
						if (p_max[n] < pop[i][n]) p_max[n] = pop[i][n];
						if (p_min[n] > pop[i][n]) p_min[n] = pop[i][n];
					}
				}
				double ms = 0;
				for (int n = 0; n < eObjNum; n++)
				{
					double lebesgue = 0;
					if (p_min[n] < pf_max[n])
					{
						double minv = pf_max[n]<p_max[n] ? pf_max[n] : p_max[n];
						double maxv = pf_min[n]>p_min[n] ? pf_min[n] : p_min[n];
						lebesgue = minv - maxv;
					}
					ms += pow(lebesgue / (pf_max[n] - pf_min[n]), 2.0);
				}
				ms_value = sqrt(ms / eObjNum);
			}

			// calculate the distance of two vectors
			double EVAL::dist_vector(vector<double> &vec1, vector<double> &vec2)
			{
				int dim = vec1.size();
				double sum = 0;
				for (int n = 0; n<dim; n++)
					sum += (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
				return sqrt(sum);
			}

			// calculate the average value of a given set
			double EVAL::Avg(vector <double> &vec)
			{
				double sum = 0.0;
				for (int i = 0; i != vec.size(); i++)
				{
					sum += vec[i];
				}
				if (vec.size()>0) return sum / vec.size();
				else return 0;
			}

			// calculate the standard deviation of a given set
			double EVAL::STDev(vector <double> &vec)
			{
				double avgv = Avg(vec);
				int size = vec.size();
				double stdv = 0;
				for (int i = 0; i < size; i++)
					stdv += pow(vec[i] - avgv, 2.0);
				if (size>1) return sqrt(stdv / (size - 1));
				else return 0;
			}


			void EVAL::Hypervolume(vector<vector<double>> pop)
			{
				// global varation 
				hv_value = 0.0;

				int pfsize = pf.size();
				vector<double> refpoint(eObjNum, -1.0e+30);
				for (int n = 0; n < eObjNum; n++)
				{
					for (int i = 0; i < pfsize; i++)
					{
						if (refpoint[n] < pf[i][n]) refpoint[n] = pf[i][n];
					}
				}
				for (int i = 0; i < refpoint.size(); i++) refpoint[i] = refpoint[i] * 1.1;

				vector<vector<double>> PopObj;
				for (int i = 0; i < pop.size(); i++)
				{
					int n = 0;
					for (; n < eObjNum; n++) if (pop[i][n] > refpoint[n]) break;
					if (n == eObjNum) PopObj.push_back(pop[i]);
				}

				if (PopObj.empty())
					return;


				sort(PopObj.begin(), PopObj.end());
				vector<pair<double, vector<vector<double>>>> S;
				S.push_back(make_pair(1, PopObj));
				for (int k = 0; k < eObjNum - 1; k++)
				{
					vector<pair<double, vector<vector<double>>>> s_;
					for (int i = 0; i < S.size(); i++)
					{
						vector<pair<double, vector<vector<double>>>> stemp;
						Hypervolume_slice(stemp, S[i].second, k, refpoint);
						for (int j = 0; j < stemp.size(); j++)
							Hypervolume_Add(make_pair(stemp[j].first*S[i].first, stemp[j].second), s_);
					}
					S = s_;
				}
				for (int i = 0; i < S.size(); i++)
					hv_value += S[i].first * fabs(S[i].second[0][eObjNum - 1] - refpoint[eObjNum - 1]);

			}


			void EVAL::Hypervolume_slice(vector<pair<double, vector<vector<double>>>> &stemp, vector<vector<double>> p1, int k, vector<double> refpoint)
			{
				vector<double> p(p1[0]);
				p1.erase(p1.begin());
				vector<vector<double>> q1;
				while (!p1.empty())
				{
					Hypervolume_Insert(p, k + 1, q1);
					Hypervolume_Add(make_pair(fabs(p[k] - p1[0][k]), q1), stemp);
					p = p1[0];
					p1.erase(p1.begin());
				}
				Hypervolume_Insert(p, k + 1, q1);
				Hypervolume_Add(make_pair(fabs(p[k] - refpoint[k]), q1), stemp);
			}


			void EVAL::Hypervolume_Insert(vector<double> p, int k, vector<vector<double>> &p1)
			{
				vector<vector<double>> q1;
				vector<double> hp;
				while (!p1.empty())
				{
					hp = p1[0];
					if (hp[k] > p[k]) break;
					q1.push_back(hp);
					p1.erase(p1.begin());
				}

				q1.push_back(p);
				while (!p1.empty())
				{
					hp = p1[0];
					for (int i = k; i < eObjNum; i++)
					{
						if (p[i] > hp[i]) { q1.push_back(hp); break; }
					}
					p1.erase(p1.begin());
				}
				p1 = q1;
			}


			void EVAL::Hypervolume_Add(pair<double, vector<vector<double>>> cell_, vector<pair<double, vector<vector<double>>>> &S)
			{
				int n = S.size();
				bool flag = false;
				for (int k = 0; k < n; k++)
				{
					if (cell_.second == S[k].second)
					{
						S[k].first = S[k].first + cell_.first;
						flag = true;
						break;
					}
				}
				if (!flag) S.push_back(cell_);
			}

		} //namespace cmea
	} //namespace mea
} //namespace az


