//	AlgD.cpp
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include "AlgC.h"
#include "Sel.h"
#include "Gen.h"

//#include "alg/AR.h"
//#include "alg/Matrix.h"
//#include "emo/GenMod.h"
//#include "alg/normal.h"
#include "assert.h"


#include <algorithm>
#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif

//#include"engine.h"
#include"memory.h"
//#pragma comment(lib,"libmat.lib") 
//#pragma comment(lib,"libmx.lib") 
//#pragma comment(lib,"libeng.lib") 
//#define SAVE_CEN 1

//!\brief	az namespace, the top namespace
namespace az
{
	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{
		//!\brief namespace of constrained evolutionary algoirhtm
		namespace cmea
		{
			const double PI = 3.141592653589793;

			// Ronsenbrock distance function
			double Ronsenbk(std::vector<double> X, unsigned int obj_num)
			{
				double sum = 1.0;
				for (unsigned int i = obj_num - 1; i < X.size() - 1; i++)
				{
					sum += 100*pow(X[i+1]-X[i]*X[i],2.0)+pow(1-X[i],2.0);
				}
				return sum;
			}

			// Bias-type distance function
			double BiasT(std::vector<double> X, unsigned int obj_num)
			{
				double sum = 1.0;
				double n = 1.0*X.size() - (double)(obj_num);
				for (unsigned int i = obj_num - 1; i < X.size(); i++)
				{
					double z = pow(X[i], n);
					sum += 5.0 - 5.0*exp(-10.0*(z - 0.5 - i / (2 * X.size()))*(z - 0.5 - i / (2 * X.size())));
				}
				return sum;
			}

			// CTP1
			void CTP1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				//double g = 1.0;
				//for (unsigned int i = 1; i < X.size(); i++)
				//	//g += X[i] * X[i];
				//	g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*exp(-1.0*F[0] / g);
				I[0] = 0.858*exp(-0.541*F[0]) - F[1];
				I[1] = 0.728*exp(-0.295*F[0]) - F[1];
			}

			// constraint for CTP2-CTP7
			double con_ctp(std::vector< double >& F, double theta, double e, double a, double b, double c, double d)
			{
				double l1 = pow((sin(theta)*(F[1] - e) + cos(theta)*F[0]), c);
				double l2 = sin(theta)*F[0] - cos(theta)*(F[1] - e);
				return a*pow(fabs(sin(b*PI*l1)), d) + l2;
			}

			// CTP2
			void CTP2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 6.0);
				//I[0] = 0.0;
			}
			
			// CTP3
			void CTP3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.1, 10.0, 1.0, 0.5);
			}
			
			// CTP4
			void CTP4(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
			}
			
			// CTP5
			void CTP5(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.1, 10.0, 2.0, 0.5);
			}
			
			// CTP6
			void CTP6(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, 0.1*PI, -2.0, 40.0, 0.5, 1.0, 2.0);
				//I[0] = 0;
			}
			
			// CTP7
			void CTP7(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, -0.05*PI, 0.0, 40.0, 5.0, 1.0, 6.0);
			}
			
			// CTP8
			void CTP8(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - F[0] / g);
				I[0] = con_ctp(F, 0.1*PI, -2.0, 40.0, 0.5, 1.0, 2.0);
				I[1] = con_ctp(F, -0.05*PI, 0.0, 40.0, 2.0, 1.0, 6.0);
			}

			// NCTP series
			// linear constraint (C2) for NCTPs
			double con_linear(std::vector< double >& F, double z)
			{
				return F[1] + 0.73*F[0] - z;
			}
			// NCTP1
			void NCTP1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);
			}

			// NCTP2
			void NCTP2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				//double g = 1.0;
				//for (unsigned int i = 1; i < X.size(); i++)
				//g += X[i] * X[i];

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);
			}

			// NCTP3
			void NCTP3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);
				I[1] = con_linear(F, 6.0);
			}

			// NCTP4
			void NCTP4(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);
			}

			// NCTP5
			void NCTP5(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
			}

			// NCTP6
			void NCTP6(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) - 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);
			}

			// NCTP7
			void NCTP7(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				//double g = 1.0;
				//for (unsigned int i = 1; i < X.size(); i++)
				//g += X[i] * X[i];

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);

			}

			// NCTP8
			void NCTP8(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);

			}

			// NCTP9
			void NCTP9(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);
				I[1] = con_linear(F, 6.0);

			}

			// NCTP10
			void NCTP10(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);

			}

			// NCTP11
			void NCTP11(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);

			}

			// NCTP12
			void NCTP12(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 0.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);

			}

			// NCTP13
			void NCTP13(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 1.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);
			}

			// NCTP14
			void NCTP14(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
				I[1] = con_linear(F, 4.0);
			}

			// NCTP15
			void NCTP15(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 3.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);
				I[1] = con_linear(F, 6.0);
			}

			// NCTP16
			void NCTP16(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 1.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.2, 10.0, 1.0, 0.5);
			}

			// NCTP17
			void NCTP17(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 1.5;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 0.75, 10.0, 1.0, 0.5);
			}

			// NCTP18
			void NCTP18(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i] - 10.0*cos(2 * PI*X[i]) + 10.0;*/

				double g;
				g = BiasT(X, F.size());

				F[0] = X[0];
				F[1] = g*(1.0 - sqrt(F[0]) / g) + 3.0;
				I[0] = con_ctp(F, -0.2*PI, 1.0, 2.0, 10.0, 1.0, 6.0);
			}

			// MW series
			// Local adjustment (LA)
			double LA1(double A, double B, double C, double D, double theta)
			{
				double t = pow(theta, C);
				return A*pow(sin(B*PI*t), D);
			}
			double LA2(double A, double B, double C, double D, double theta)
			{
				double t = pow(theta, C);
				return A*pow(sin(B*t), D);
			}
			double LA3(double A, double B, double C, double D, double theta)
			{
				double t = pow(theta, C);
				return A*pow(cos(B*t), D);
			}

			// distance functions
			// g1
			double g1(std::vector<double> X, unsigned int obj_num)
			{
				double sum = 1.0;
				double n = 1.0*X.size()-(double)(obj_num);
				for (unsigned int i = obj_num - 1; i < X.size(); i++)
				{	
					double z = pow(X[i],n);
					sum += 1 - exp(-10.0*(z - 0.5 - i / (2 * X.size()))*(z - 0.5 - i / (2 * X.size())));
				}
				return sum;
			}

			// g2
			double g2(std::vector<double> X, unsigned int obj_num)
			{
				double sum = 1.0;
				double n = 1.0*X.size();
				for (unsigned int i = obj_num - 1; i < X.size(); i++)
				{
					double z = 1 - exp(-10.0*(X[i] - i / n)*(X[i] - i / n));
					sum += (0.1/(n))*z*z + 1.5 - 1.5*cos(2 * PI*z);
				}
				return sum;
			}

			// g3
			double g3(std::vector<double> X, unsigned int obj_num)
			{
				double sum = 1.0;
				double n = 1.0*X.size();
				for (unsigned int i = obj_num - 1; i < X.size(); i++)
				{
					sum += 2.0 * pow(X[i] + (X[i - 1] - 0.5)*(X[i - 1] - 0.5) - 1.0, 2.0);
				}
				return sum;
			}

			// MW1:
			void CMP1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g1(X, F.size());

				F[0] = X[0];
				F[1] = g*(1 - 0.85*F[0] / g);
				I[0] = F[0] + F[1] - 1 - LA1(0.5, 2.0, 1.0, 8.0, sqrt(2.0)*F[1] - sqrt(2.0)*F[0]);

			}

			// MW2
			void CMP2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g2(X, F.size());

				F[0] = X[0];
				F[1] = g*(1 - F[0] / g);
				I[0] = (F[0] + F[1] - 1 - LA1(0.5, 3.0, 1.0, 8.0, sqrt(2.0)*F[1] - sqrt(2.0)*F[0]));
			}

			// MW3
			void CMP3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g3(X, F.size());

				F[0] = X[0];
				F[1] = g*(1 - F[0] / g);
				I[0] = F[0] + F[1] - 1.05 - LA1(0.45, 0.75, 1.0, 6.0, sqrt(2.0)*F[1] - sqrt(2.0)*F[0]);
				I[1] = 0.85 - F[0] - F[1] + LA1(0.3, 0.75, 1.0, 2.0, sqrt(2.0)*F[1] - sqrt(2.0)*F[0]);
			}

			// MW4
			void CMP4(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 2; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g1(X, F.size());

				F[0] = g* (1-X[0]) * (1-X[1]);
				F[1] = g* (1-X[0]) * (X[1]);
				F[2] = g*(X[0]);
				I[0] = F[0] + F[1] + F[2] - (1.0 + LA1(0.4, 2.5, 1.0, 8.0, F[2] - F[1] - F[0]));
				//I[0] = 0;
			}

			// MW5
			void CMP5(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g1(X, F.size());

				F[0] = g*X[0];
				F[1] = g*sqrt(1.0-pow(F[0]/g,2.0));
				I[0] = F[0] * F[0] + F[1] * F[1] - pow(1.7 - LA2(0.2, 2.0, 1.0, 1.0, atan(F[1] / F[0])), 2.0);
				double t = 0.5*PI - 2 * fabs(atan(F[1] / F[0]) - 0.25*PI);
				I[1] = pow(1 + LA2(0.5, 6.0, 3.0, 1.0, t), 2.0) - F[0] * F[0] - F[1] * F[1];
				I[2] = pow(1 - LA2(0.45, 6.0, 3.0, 1.0, t), 2.0) - F[0] * F[0] - F[1] * F[1];
			}
			
			// MW6
			void CMP6(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g2(X, F.size());

				F[0] = g*X[0];
				F[1] = g*sqrt(1.1*1.1 - pow(F[0] / g, 2.0));
				I[0] = F[0] * F[0] / pow(1.0 + LA3(0.15, 6.0, 4.0, 10.0, atan(F[1] / F[0])), 2.0) + F[1] * F[1] / pow(1.0 + LA3(0.75, 6.0, 4.0, 10.0, atan(F[1] / F[0])), 2.0) - 1;
			}

			// MW7
			void CMP7(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i];*/

				double g;
				g = g3(X, F.size());

				F[0] = g*X[0];
				F[1] = g*sqrt(1 - pow(F[0] / g, 2.0));
				I[0] = F[0] * F[0] + F[1] * F[1] - pow(1.2 + fabs(LA2(0.4, 4.0, 1.0, 16.0, atan(F[1] / F[0]))), 2.0);
				I[1] = pow(1.15 - LA2(0.2, 4.0, 1.0, 8.0, atan(F[1] / F[0])), 2.0) - F[0] * F[0] - F[1] * F[1];
			}

			// MW8
			void CMP8(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 2; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g2(X, F.size());

				F[0] = g*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
				F[1] = g*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
				F[2] = g*sin(0.5*PI*X[0]);
				I[0] = F[0] * F[0] + F[1] * F[1] + F[2] * F[2] - (1.25 - LA2(0.5, 6.0, 1.0, 2.0, asin(F[2])))*(1.25 - LA2(0.5, 6.0, 1.0, 2.0, asin(F[2])));
			}

			// MW9
			void CMP9(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i];*/

				double g;
				g = g1(X, F.size());
				
				F[0] = g*X[0];
				F[1] = g*(1.0 - pow(F[0] / g, 0.6));

				double T1 = (1 - 0.64*F[0] * F[0] - F[1])*(1 - 0.36*F[0] * F[0] - F[1]);
				double T2 = (1.35*1.35 - (F[0] + 0.35)*(F[0] + 0.35) - F[1])*(1.15*1.15 - (F[0] + 0.15)*(F[0] + 0.15) - F[1]);
				I[0] = T1 <= T2 ? T1 : T2;
			}

			// MW10
			void CMP10(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g2(X, F.size());

				F[0] = g*pow(X[0], X.size());
				F[1] = g*(1.0 - pow(F[0] / g, 2.0));
				I[0] = -1.0*(2.0 - 4.0*F[0] * F[0] - F[1])*(2.0 - 8.0*F[0] * F[0] - F[1]);
				I[1] = (2.0 - 2.0*F[0] * F[0] - F[1])*(2.0 - 16.0*F[0] * F[0] - F[1]);
				I[2] = (1.0 - F[0] * F[0] - F[1])*(1.2 - 1.2*F[0] * F[0] - F[1]);
			}

			// MW11
			void CMP11(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g3(X, F.size());

				F[0] = g*X[0];
				F[1] = g*sqrt(2.0 - pow(F[0] / g, 2.0));
				I[0] = -1.0*(3.0 - F[0] * F[0] - F[1])*(3.0 - 2.0*F[0] * F[0] - F[1]);
				I[1] = (3.0 - 0.625*F[0] * F[0] - F[1])*(3.0 - 7.0*F[0] * F[0] - F[1]);
				I[2] = -1.0*(1.62 - 0.18*F[0] * F[0] - F[1])*(1.125 - 0.125*F[0] * F[0] - F[1]);
				I[3] = (2.07 - 0.23*F[0] * F[0] - F[1])*(0.63 - 0.07*F[0] * F[0] - F[1]);
			}

			// MW12
			void CMP12(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
				g += X[i] * X[i];*/

				double g;
				g = g1(X, F.size());

				F[0] = g*X[0];
				F[1] = g*(0.85 - 0.8*(F[0] / g) - 0.08*fabs(sin(3.2*PI*(F[0]/g))));
				I[0] = -1.0*(1 - 0.625*F[0] - F[1] + 0.08*sin(2 * PI*(F[1] - F[0] / 1.6)))*(1.4 - 0.875*F[0] - F[1] + 0.08*sin(2 * PI*(F[1]/1.4 - F[0] / 1.6)));
				I[1] = (1 - 0.8*F[0] - F[1] + 0.08*sin(2 * PI*(F[1] - F[0] / 1.5)))*(1.8 - 1.125*F[0] - F[1] + 0.08*sin(2 * PI*(F[1]/1.8 - F[0] / 1.6)));
			}

			// MW13
			void CMP13(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 1; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g2(X, F.size());

				F[0] = g*X[0];
				F[1] = g*(5.0 - exp(F[0] / g) - fabs(0.5*sin(3 * PI*F[0] / g)));
				I[0] = -1.0*(5.0 - (1 + F[0] + 0.5*F[0] * F[0]) - 0.5*sin(3 * PI*F[0]) - F[1])*(5.0 - (1 + 0.7*F[0]) - 0.5*sin(3 * PI*F[0]) - F[1]);
				I[1] = (5.0 - exp(F[0]) - 0.5*sin(3 * PI*F[0]) - F[1])*(5.0 - (1 + 0.4*F[0]) - 0.5*sin(3 * PI*F[0]) - F[1]);
			}

			// MW14
			void CMP14(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
			{
				/*double g = 1.0;
				for (unsigned int i = 2; i < X.size(); i++)
					g += X[i] * X[i];*/

				double g;
				g = g3(X, F.size());

				F[0] = X[0];
				F[1] = X[1];
				F[2] = 0.5*g*((6 - exp(F[0]) - LA1(1.5, 1.1, 2.0, 1.0, F[0])) + (6 - exp(F[1]) - LA1(1.5, 1.1, 2.0, 1.0, F[1])));
				//I[0] = (6 - exp(F[0]) - LA1(1.5, 1.1, 2.0, 1.0, F[0])) + (6 - exp(F[1]) - LA1(1.5, 1.1, 2.0, 1.0, F[1]));
				I[0] = F[2]-0.5*(6.1 - (1 + F[0] + 0.5*F[0] * F[0]) - LA1(1.5, 1.1, 2.0, 1.0, F[0]) + 6.1 - (1 + F[1] + 0.5*F[1] * F[1]) - LA1(1.5, 1.1, 2.0, 1.0, F[1]));
				//I[0] = 0.0;
			}


			void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name, unsigned int mObj)
			{
				if (  name == std::string("CTP1") || name == std::string("CTP2") || name == std::string("CTP3") || name == std::string("CTP4")
						|| name == std::string("CTP5") || name == std::string("CTP6") || name == std::string("CTP7") || name == std::string("CTP8"))
				{
					low[0] = 0.0; upp[0] = 1.0;
					for (unsigned int i = 1; i<(unsigned int)(low.size()); i++){ low[i] = -1.0; upp[i] = 1.0; }
				}
				else if (  name == std::string("CMP1") || name == std::string("CMP2") || name == std::string("CMP3") || name == std::string("CMP4") 
						|| name == std::string("CMP5") || name == std::string("CMP7") || name == std::string("CMP8") || name == std::string("CMP9")
						|| name == std::string("CMP10")	|| name == std::string("CMP12"))
				{
					for (unsigned int i = 0; i < (unsigned int)(low.size()); i++){ low[i] = 0.0; upp[i] = 1.0; }
				}
				else if (name == std::string("CMP6"))
				{
					for (unsigned int i = 0; i < (unsigned int)(low.size()); i++){ low[i] = 0.0; upp[i] = 1.1; }
				}
				else if (name == std::string("CMP11"))
				{
					for (unsigned int i = 0; i < (unsigned int)(low.size()); i++){ low[i] = 0.0; upp[i] = sqrt(2); }
				}
				else if (name == std::string("CMP13") || name == std::string("CMP14"))
				{
					for (unsigned int i = 0; i < (unsigned int)(low.size()); i++){ low[i] = 0.0; upp[i] = 1.5; }
				}
				else if (name == std::string("NCTP1") || name == std::string("NCTP2") || name == std::string("NCTP3") || name == std::string("NCTP4")
					|| name == std::string("NCTP5") || name == std::string("NCTP6") || name == std::string("NCTP7") || name == std::string("NCTP8")
					|| name == std::string("NCTP9") || name == std::string("NCTP10") || name == std::string("NCTP11") || name == std::string("NCTP12")
					|| name == std::string("NCTP13") || name == std::string("NCTP14") || name == std::string("NCTP15") || name == std::string("NCTP16")
					|| name == std::string("NCTP17") || name == std::string("NCTP18"))
				{
					low[0] = 0.0; upp[0] = 5.0;
					for (unsigned int i = 1; i<(unsigned int)(low.size()); i++){ low[i] = 0.0; upp[i] = 1.0; }
				}
			}

			CMOO::CMOO(
				std::string&	optimizer,
				unsigned int	popsize,
				unsigned int	stepmax,
				CParameter&	par)
				:mPop(par), ex_pop(par)
			{
				mOptimizer = optimizer;
				mPopSize = popsize;
				mMaxStep = stepmax;
				pPar = &par;

				
				if (pPar->Problem() == std::string("CTP1"))
				{
					P().Evaluator(CTP1);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("CTP2"))
				{
					P().Evaluator(CTP2);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP3"))
				{
					P().Evaluator(CTP3);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP4"))
				{
					P().Evaluator(CTP4);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP5"))
				{
					P().Evaluator(CTP5);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP6"))
				{
					P().Evaluator(CTP6);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP7"))
				{
					P().Evaluator(CTP7);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CTP8"))
				{
					P().Evaluator(CTP8);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP1"))
				{
					P().Evaluator(NCTP1);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP2"))
				{
					P().Evaluator(NCTP2);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP3"))
				{
					P().Evaluator(NCTP3);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP4"))
				{
					P().Evaluator(NCTP4);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP5"))
				{
					P().Evaluator(NCTP5);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP6"))
				{
					P().Evaluator(NCTP6);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP7"))
				{
					P().Evaluator(NCTP7);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP8"))
				{
					P().Evaluator(NCTP8);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP9"))
				{
					P().Evaluator(NCTP9);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP10"))
				{
					P().Evaluator(NCTP10);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP11"))
				{
					P().Evaluator(NCTP11);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP12"))
				{
					P().Evaluator(NCTP12);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP13"))
				{
					P().Evaluator(NCTP13);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP14"))
				{
					P().Evaluator(NCTP14);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP15"))
				{
					P().Evaluator(NCTP15);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("NCTP16"))
				{
					P().Evaluator(NCTP16);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP17"))
				{
					P().Evaluator(NCTP17);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("NCTP18"))
				{
					P().Evaluator(NCTP18);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP1"))
				{
					P().Evaluator(CMP1);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP2"))
				{
					P().Evaluator(CMP2);
					P().FSize(2); P().ESize(1); P().ISize(0);
				}
				else if (pPar->Problem() == std::string("CMP3"))
				{
					P().Evaluator(CMP3);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("CMP4"))
				{
					P().Evaluator(CMP4);
					P().FSize(3); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP5"))
				{
					P().Evaluator(CMP5);
					P().FSize(2); P().ESize(0); P().ISize(3);
				}
				else if (pPar->Problem() == std::string("CMP6"))
				{
					P().Evaluator(CMP6);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP7"))
				{
					P().Evaluator(CMP7);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("CMP8"))
				{
					P().Evaluator(CMP8);
					P().FSize(3); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP9"))
				{
					P().Evaluator(CMP9);
					P().FSize(2); P().ESize(0); P().ISize(1);
				}
				else if (pPar->Problem() == std::string("CMP10"))
				{
					P().Evaluator(CMP10);
					P().FSize(2); P().ESize(0); P().ISize(3);
				}
				else if (pPar->Problem() == std::string("CMP11"))
				{
					P().Evaluator(CMP11);
					P().FSize(2); P().ESize(0); P().ISize(4);
				}
				else if (pPar->Problem() == std::string("CMP12"))
				{
					P().Evaluator(CMP12);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("CMP13"))
				{
					P().Evaluator(CMP13);
					P().FSize(2); P().ESize(0); P().ISize(2);
				}
				else if (pPar->Problem() == std::string("CMP14"))
				{
					P().Evaluator(CMP14);
					P().FSize(3); P().ESize(0); P().ISize(1);
				}

				Reset();
			}

			void CMOO::Reset()
			{
				mStep = 0;
				mEvas = 0;
				mCP = 2.0;

				ex_pop.Clear();		// clear the external population

				DFRange(P().XLow(), P().XUpp(), P().Problem(), P().FSize());	// set the original search space

				az::rnd::seed((long)time(NULL));
			}

			// main evolution step
			unsigned int CMOO::Step()
			{
				// initialization of population
				if (CurStep() == 0) 
				{
					Init(); Population().Evaluate(); mEvas += mPopSize;
				}

				// evolution
				CPopulationMO popnew(P());
				// generate new solutions
				Generate(popnew, mPopSize);

				//check the range of new solutions
				Check(popnew);
				//remove these repeated solutions
				Check(popnew, Population());

				//evaluate new solutions
				popnew.Evaluate(); mEvas += popnew.Size();
				
				//environmental select
				Population().Combine(popnew);
				Select(Population(), mPopSize);

				mStep++;
				return mStep;
			}

			// initialization
			void CMOO::Init()
			{
				unsigned int i, j;
				Population().Resize(mPopSize);
				for (i = 0; i < mPopSize; i++) for (j = 0; j < P().XSize(); j++)
				{
					Population()[i][j] = rnd::rand(P().XLow(j), P().XUpp(j));
					Population()[i].mID = i;
				}
			}

			// check the range of each individual
			void CMOO::Check(CPopulationMO& pop)
			{
				for (unsigned int i = 0; i<pop.Size(); i++) pop[i].Check();
			}

			// delete duplicate individuals
			void CMOO::Check(CPopulationMO& popnew, CPopulationMO& pop)
			{
				CPopulationMO tmp(P());
				for (unsigned int i = 0; i<popnew.Size(); i++) if (!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
				popnew.Clear();
				popnew.Combine(tmp);
			}

			// generate new trial solutions by any offspring generator
			CPopulationMO& CMOO::Generate(CPopulationMO& popnew, unsigned int size)
			{
				if (mOptimizer == std::string("ToR"))
				{
					popnew.P().ETA_PM() = 20.0;
					popnew.P().ETA_SBX() = 20.0;
					popnew.P().Pc() = 0.9;
					popnew.P().Pm() = 1.0 / popnew.P().XSize();
					az::mea::gen::XSBX gen;
					gen.Generate(size, popnew, Population());

					//az::mea::gen::XDE gen;
					//gen.Set(0.4, 0.8);
					//gen.Generate(size, popnew, Population());
				}
				return popnew;
			}

			//environmental selection
			CPopulationMO& CMOO::Select(CPopulationMO& pop, unsigned int size)
			{
				if (mOptimizer == std::string("ToR"))
				{
					Population().P().EPSN_CHT() = (1.0*mStep) / mMaxStep;
					az::mea::sel::SDouCrowd sel;
					sel.Select(pop, size);
				}
				
				// get the number of feasible solutions in pop
				pop.FN() = 0;
				for (unsigned int i = 0; i < pop.Size(); i++) if (pop[i].IsFeasible()) pop.FN() = pop.FN() + 1;

				return pop;
			}

			// Plot function
			void ePlot(Engine*& ep, std::string& problem, unsigned int obj)
			{
				if (problem == std::string("CTP1"))
				{
					engEvalString(ep, "frontx=0:0.01:1/3; fronty1=exp(-1.0*frontx); fronty2=0.858*exp(-0.541*(frontx+1/3)); fronty3=0.728*exp(-0.295*(frontx+2/3));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "plot(frontx,fronty1,'Color',[0 0.45 0.74],'Linewidth',2); plot(frontx+1/3,fronty2,'Color',[0.64 0.08 0.18],'Linewidth',1.5); plot(frontx+2/3,fronty3,'Color',[0.64 0.08 0.18],'Linewidth',1.5)");
				}
				else if (problem == std::string("CTP2"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.2*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^6+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,1,0,1.5]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("CTP3"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.1*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,1,0,1.5]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("CTP4"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,500),linspace(0,1,500));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h0=ezplot('sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,1]); set(h0,'LineWidth',0.5,'LineColor',[0.75 0.75 0.75],'LineStyle',':')");
					engEvalString(ep, "h1=ezplot('0.75*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,1,0,2]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("CTP5"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.1*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^2)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,1,0,1.5]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");

				}
				else if (problem == std::string("CTP6"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('40*abs(sin(0.5*pi*(sin(0.1*pi)*(y+2)+cos(0.1*pi)*x).^1)).^2+sin(0.1*pi)*x-cos(0.1*pi)*(y+2)',[0,1,0,20]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("CTP7"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('40*abs(sin(5*pi*(sin(-0.05*pi)*(y)+cos(-0.05*pi)*x).^1)).^6+sin(-0.05*pi)*x-cos(-0.05*pi)*(y)',[0,1,0,2]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("CTP8"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h0=ezplot('40*abs(sin(0.5*pi*(sin(0.1*pi)*(y+2)+cos(0.1*pi)*x).^1)).^2+sin(0.1*pi)*x-cos(0.1*pi)*(y+2)',[0,1,0,20]); set(h0,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]);");
					engEvalString(ep, "h1=ezplot('40*abs(sin(2*pi*(sin(-0.05*pi)*(y)+cos(-0.05*pi)*x).^1)).^6+sin(-0.05*pi)*x-cos(-0.05*pi)*(y)',[0,1,0,20]); h2=ezplot('1-x-y',[0,1]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP1") || problem == std::string("NCTP4"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.2*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y-1.5',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP2") || problem == std::string("NCTP5"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.75*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y-1.5',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP3") || problem == std::string("NCTP6"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('2.0*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^6+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y-1.5',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP7") || problem == std::string("NCTP10"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.2*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP8") || problem == std::string("NCTP11"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.75*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP9") || problem == std::string("NCTP12"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('2.0*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^6+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('1-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP13") || problem == std::string("NCTP16"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.2*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('2-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP14") || problem == std::string("NCTP17"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('0.75*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^0.5+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('2.5-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				else if (problem == std::string("NCTP15") || problem == std::string("NCTP18"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,5,5000),linspace(0,5,5000));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "h1=ezplot('2.0*abs(sin(10*pi*(sin(-0.2*pi)*(y-1)+cos(-0.2*pi)*x).^1)).^6+sin(-0.2*pi)*x-cos(-0.2*pi)*(y-1)',[0,5,-3,4]); h2=ezplot('4-sqrt(x)-y',[0,5,-3,4]); set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74])");
				}
				
				else if (problem == std::string("CMP1"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('1-x-y+0.5*sin(2.0*pi*sqrt(2.0)*(y-x)).^8',[0,1,0,2]); h2=ezplot('1-0.85*x-y',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP1');");
				}
				else if (problem == std::string("CMP2"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('1-x-y+0.5*sin(3.0*pi*sqrt(2.0)*(y-x)).^8',[0,1,0,2]); h2=ezplot('1-x-y',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP2')");
				}
				else if (problem == std::string("CMP3"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('1.05-x-y+0.45*sin(0.75*pi*sqrt(2.0)*(y-x)).^6',[0,1,0,2]); h2=ezplot('0.85-x-y+0.3*sin(0.75*pi*sqrt(2.0)*(y-x)).^2',[0,1,0,2]); h3=ezplot('1-x-y',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h3,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP3')");
				}
				else if (problem == std::string("CMP4"))
				{
					engEvalString(ep, "[t1,t2]=meshgrid(linspace(0,1,50),linspace(0,1,50)); frontx=t1.*t2; fronty=t1.*(1-t2); frontz=(1-t1);");
					engEvalString(ep, "hold on");
					engEvalString(ep, "surf(frontx,fronty,frontz,'FaceAlpha',0.6,'FaceColor',[0 0.45 0.74],'EdgeColor','none');");
					engEvalString(ep, "cx=frontx.*(1.01+0.4*sin(2.5*pi*(frontz-fronty-frontx)).^8); cy=fronty.*(1.01+0.4*sin(2.5*pi*(frontz-fronty-frontx)).^8); cz=frontz.*(1.01+0.4*sin(2.5*pi*(frontz-fronty-frontx)).^8);");
					engEvalString(ep, "surf(cx,cy,cz,'FaceAlpha',0.3,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')");
				}
				else if (problem == std::string("CMP5"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(1.7-0.2*sin(2*atan(y./x))).^2-x.^2-y^2',[0,2,0,2]); h2=ezplot('max((1.0+0.5*sin(6*(pi/2-2*abs(atan(y./x)-pi/4)).^3)).^2-x.^2-y^2,1-x.^2-y.^2)',[0,2]); h3=ezplot('max((1.0-0.45*sin(6*(pi/2-2*abs(atan(y./x)-pi/4)).^3)).^2-x.^2-y^2,1-x.^2-y.^2)',[0,2]);h4=ezplot('1-x.^2-y^2',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h3,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h4,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP5')");
				}
				else if (problem == std::string("CMP6"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(x./(1+0.15*cos(6*atan(y./x).^4).^10)).^2+(y./(1+0.75*cos(6*atan(y./x).^4).^10)).^2-1',[0,2,0,2]); h2=ezplot('1.1*1.1-x.^2-y^2',[0,1.1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP6')");
				}
				else if (problem == std::string("CMP7"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(1.2+0.4*abs(sin(4*atan(y./x))).^16).^2-x.^2-y^2',[0,2,0,2]); h2=ezplot('(1.15-0.2*sin(4*atan(y./x)).^8).^2-x.^2-y^2',[0,2,0,2]); h3=ezplot('1-x.^2-y^2',[0,1.1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h3,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP7')");
				}
				else if (problem == std::string("CMP8"))
				{
					engEvalString(ep, "[t1,t2]=meshgrid(linspace(0,1,100),linspace(0,1,100)); frontx=cos(0.5*pi*t1).*cos(0.5*pi*t2); fronty=cos(0.5*pi*t1).*sin(0.5*pi*t2); frontz=sin(0.5*pi*t1);");
					engEvalString(ep, "hold on");
					engEvalString(ep, "surf(frontx,fronty,frontz,'FaceAlpha',0.6,'FaceColor',[0 0.45 0.74],'EdgeColor','none');");
					engEvalString(ep, "cx=frontx.*(1.25-0.5*sin(6*asin(frontz)).^2); cy=fronty.*(1.25-0.5*sin(6*asin(frontz)).^2); cz=frontz.*(1.25-0.5*sin(6*asin(frontz)).^2);");
					engEvalString(ep, "surf(cx,cy,cz,'FaceAlpha',0.4,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')");
				}
				else if (problem == std::string("CMP9"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(1-0.64*x.^2-y).*(1-0.36*x.^2-y).*(1.35^2-(x+0.35).^2-y).*(1.15^2-(x+0.15).^2-y)',[0,1.8,0,1.8]); h2=ezplot('1-x.^0.6-y',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP9')");
				}
				else if (problem == std::string("CMP10"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(2-4*x.^2-y).*(2-8*x.^2-y).*(2-2*x.^2-y).*(2-16*x.^2-y)*(1-x.^2-y)*(1.2-1.2*x.^2-y)',[0,1,0,2]); h2=ezplot('1-x.^2-y',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP9')");
				}
				else if (problem == std::string("CMP11"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(3-x.^2-y).*(3-2*x.^2-y).*(3-0.625*x.^2-y).*(3-7*x.^2-y)*(1.62-0.18*x.^2-y)*(1.125-0.125*x.^2-y)*(2.07-0.23*x.^2-y)*(0.63-0.07*x.^2-y)',[0,3,0,3]); h2=ezplot('2-x.^2-y.^2',[0,2]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP10')");
				}
				else if (problem == std::string("CMP12"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,100),linspace(0,1,100));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(1-0.8*x-y+0.08*sin(2*pi*(y-x/1.5))).*(1-0.625*x-y+0.08*sin(2*pi*(y-x/1.6))).*(1.4-0.875*x-y+0.08*sin(2*pi*(y/1.4-x/1.6))).*(1.8-1.125*x-y+0.08*sin(2*pi*(y/1.8-x/1.6)))',[0,1.8,0,2]); h2=ezplot('0.85-0.8*x-y-0.08*abs(sin(3.2*pi*x))',[0,1]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP9')");
				}
				else if (problem == std::string("CMP13"))
				{
					engEvalString(ep, "[x, y]= meshgrid(linspace(0,1,1000),linspace(0,1,1000));");
					engEvalString(ep, "hold on;");
					engEvalString(ep, "h1=ezplot('(5-1-x-0.5*x.^2-0.5*sin(3*pi*x)-y).*(5-1-0.7*x-0.5*sin(3*pi*x)-y).*(5-1-0.4*x-0.5*sin(3*pi*x)-y).*(5-exp(x)-0.5*sin(3*pi*x)-y)',[0,2,0,5]); h2=ezplot('5-exp(x)-abs(0.5*sin(3*pi*x))-y',[0,2,0,5]);");
					engEvalString(ep, "set(h1,'LineWidth',1.5,'LineColor',[0.64 0.08 0.18]); set(h2,'LineWidth',1.5,'LineColor',[0 0.45 0.74]); title('CMP11')");
				}
				else if (problem == std::string("CMP14"))
				{
					engEvalString(ep, "[t1,t2]=meshgrid(linspace(0,1.5,150),linspace(0,1.5,150)); frontx=t1; fronty=t2; frontz=0.5*(6-exp(frontx)-1.5*sin(1.1*pi*frontx.^2)+6-exp(fronty)-1.5*sin(1.1*pi*fronty.^2));");
					engEvalString(ep, "hold on");
					engEvalString(ep, "surf(frontx,fronty,frontz,'FaceAlpha',0.6,'FaceColor',[0 0.45 0.74],'EdgeColor','none');");
					engEvalString(ep, "cx=t1; cy=t2; cz=0.5*(6.2-(1+cx+0.5*cx.^2)-1.5*sin(1.1*pi*cx.^2)+6.2-(1+cy+0.5*cy.^2)-1.5*sin(1.1*pi*cy.^2));");
					engEvalString(ep, "surf(cx,cy,cz,'FaceAlpha',0.4,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')");
				}

			}

			// Print Copyright Information
			void PrintCopyright()
			{
				std::cout << " ====================================================================================" << std::endl;
				std::cout << " |This is a software framework for constrained multi-objective optimization. It is  |" << std::endl;
				std::cout << " |modifeid from the dynamic framework created by Aimin Zhou <amzhou@cs.ecnu.edu.cn>.|" << std::endl;
				std::cout << " |We sincerely thank him for providing us with the source code.                     |" << std::endl;
				std::cout << " |                                                                                  |" << std::endl;
				std::cout << " |* This code is free for research pruposes ONLY! *                                 |" << std::endl;
				std::cout << " |  Copyright: Copyright (c) 2017-2018 CITAA Group                                  |" << std::endl;
				std::cout << " |  Creation date: 2017-5-1                                                         |" << std::endl;
				std::cout << " |  Author: Zhongwei Ma <mzw_cemo@csu.edu.cn>                                       |" << std::endl;
				std::cout << " |  Version: 1.0                                                                    |" << std::endl;
				std::cout << " |                                                                                  |" << std::endl;
				std::cout << " |(If you have any questions, please contact me. This is very important for our impr|" << std::endl;
				std::cout << " |-ovement. Finally, good luck to you!)                                             |" << std::endl;
				std::cout << " ====================================================================================" << std::endl << std::endl;
				std::cout << " The experiment starts..." << std::endl << std::endl;
			}

			// Print parameter settings
			void ParaSetting(std::string method, unsigned int runs, unsigned int maxgen, unsigned int popsize, unsigned int varsize)
			{
				std::cout << "  Algorithm: " << method.c_str() << "	Total Runs: " << runs << "		Max Gen: " << maxgen << std::endl;
				std::cout << "  Pop Size: " << popsize << "		Num of Var: " << varsize << std::endl << std::endl;
			}

			// Show the Rate of Progress
			void ShowProgress(unsigned int run, unsigned int gen, unsigned int ic, unsigned int maxgen)
			{
				if ((gen + 1) == ic)
				{
					std::cout << "¨€ " << (int)((gen + 1)*4.0 / ic) << "%";
				}
				else if ((gen + 1) % ic == 0 && (gen + 1) != ic)
				{
					if (((gen + 1)*4.0 / ic) <= 12)
						std::cout << "\b\b\b¨€ " << (int)((gen + 1)*4.0 / ic) << "%";
					else std::cout << "\b\b\b\b¨€ " << (int)((gen + 1)*4.0 / ic) << "%";
				}
				if ((gen + 1) == maxgen)
					std::cout << "\r\t\t\t\t\t\t\t\t\t\r  Run " << ++run << ": ";
			}


		} //namespace cmea
	} //namespace mea
} //namespace az
