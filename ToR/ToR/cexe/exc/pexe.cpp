#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "AlgC.h"
#include "Metric.h"
#include "Parameter.h"
#include "Config.h"
#include "stdlib.h"
#include "Random.h"
#include <random>

#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif

int main(int argc, char* argv[])
{
	int pro_index = 0;
	az::mea::cmea::PrintCopyright();
	while (pro_index < 14)
	{
		std::string method, problem, fname, path;
		unsigned int strategy, runs, max_gen, popsize, dimension, generation, ir, i, j, nt, ic;
		char savefilename[1024];
		bool DebugModel;
		clock_t start, finish;

		//Parameter setting 
		std::string CTP[] = { "CTP1", "CTP2", "CTP3", "CTP4", "CTP5", "CTP6", "CTP7", "CTP8" };
		std::string NCTP[] = { "NCTP1", "NCTP2", "NCTP3", "NCTP4", "NCTP5", "NCTP6", "NCTP7", "NCTP8", "NCTP9", "NCTP10", "NCTP11", "NCTP12", "NCTP13", "NCTP14", "NCTP15", "NCTP16", "NCTP17", "NCTP18" };
		// MW test problems published recently on TEVC by the authors:
		std::string CMP[] = { "CMP1", "CMP2", "CMP3", "CMP4", "CMP5", "CMP6", "CMP7", "CMP8", "CMP9", "CMP10", "CMP11", "CMP12", "CMP13", "CMP14" };

		std::string arraymethod[] = { "ToR" };

		problem = CMP[pro_index];

		method = arraymethod[0];

		runs = 5;						//Number of runs
		popsize = 100;
		dimension = 15;
		max_gen = 500;

		// true or flase
		// if true, debug model will call matlab to plot the poppulation distribution
		// (it requires matlab version of 2014a and system path should be correctly configured!)
		DebugModel = false;
		nt = 20;

		az::mea::CParameter mPar;
		mPar.TolF() = 0.0;				//TOLERANCEF
		mPar.TolX() = 1.0E-5;				//TOLERANCEX
		mPar.TolC() = 1.0E-5;				//TOLERANCEC
		mPar.EPSN_CHT() = 0.0;
		mPar.XSize(dimension);
		mPar.Problem(problem);
		mPar.XCoding() = false;

		//unsigned int ipf;
		//std::vector< std::vector<double> > PFs, PSs;
		//std::ofstream f1, f2; 

		std::ofstream pf;

		az::mea::cmea::CMOO* pEA = new az::mea::cmea::CMOO(method, popsize, max_gen, mPar);
		az::mea::cmea::EVAL* pEval = new az::mea::cmea::EVAL(runs, max_gen, popsize, pEA->P().FSize(), problem);
		unsigned int xdim = (dimension < 3) ? dimension : 3;
		static bool flag = false;

		if (pro_index == 0) { az::mea::cmea::ParaSetting(method, runs, max_gen, popsize, dimension); }
		std::cout << " ---------------------------------------- " << problem.c_str() << " ---------------------------------------- " << std::endl << std::endl << "  Run 1: ";
		start = clock();

		for (ir = 1; ir <= runs; ir++)
		{

			/*ipf = 0;
			PFs.resize(int(popsize * max_gen * 1.0 / nt));
			PSs.resize(int(popsize * max_gen * 1.0 / nt));*/

			ic = (int)(max_gen*1.0 / 25);
			
			pEA->Reset();

			while (!pEA->IsTerminate())
			{
				Engine *ep;
				generation = pEA->Step();

				if ((generation + 1) % nt == 0 && DebugModel)
				{
					if (pEA->Population().P().FSize()==2)
					{
						double f1[100];
						double f2[100];
						double mNt = nt;
						mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
						mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
						mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
						for(int n=0; n<pEA->Population().Size(); n++)
						{ 
							
							/*f1[n] = tmp_pop[n].F(0);
							f2[n] = tmp_pop[n].F(1);*/
							f1[n] = pEA->Population()[n].F(0);
							f2[n] = pEA->Population()[n].F(1);
						}
						if(  (ep=engOpen(NULL)) )
						{
							memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
							memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
							memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
							engPutVariable(ep , "TX", T) ;
							engPutVariable(ep , "TY", T2) ;
							engPutVariable(ep , "NT", M) ;
							if(flag==false)
							{
								flag =true ;
								engEvalString( ep ,"h=plot(TX,TY,'ko');box on");
								az::mea::cmea::ePlot(ep, problem, pEA->Population().P().FSize());
							}
							else
							{
								engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY   )");
							}
						}
					}
					// 3-OBJECTIVE
					if (pEA->Population().P().FSize() == 3)
					{
						double f1[200];
						double f2[200];
						double f3[200];
						double mNt = nt;
						mxArray *T = mxCreateDoubleMatrix(1, pEA->Population().Size(), mxREAL);
						mxArray *T2 = mxCreateDoubleMatrix(1, pEA->Population().Size(), mxREAL);
						mxArray *T3 = mxCreateDoubleMatrix(1, pEA->Population().Size(), mxREAL);
						mxArray *M = mxCreateDoubleMatrix(1, 1, mxREAL);
						int c = 0;
						for (int n = 0; n < pEA->Population().Size(); n++)if (pEA->Population()[n].IsFeasible())
						{
							/*f1[n] = tmp_pop[n].F(0);
							f2[n] = tmp_pop[n].F(1);
							f3[n] = tmp_pop[n].F(2);*/
							f1[c] = pEA->Population()[n].F(0);
							f2[c] = pEA->Population()[n].F(1);
							f3[c] = pEA->Population()[n].F(2);
							c++;
						}
						if ((ep = engOpen(NULL)))
						{
							memcpy((char*)mxGetPr(T), (char*)f1, pEA->Population().Size()*sizeof(double));
							memcpy((char*)mxGetPr(T2), (char*)f2, pEA->Population().Size()*sizeof(double));
							memcpy((char*)mxGetPr(T3), (char*)f3, pEA->Population().Size()*sizeof(double));
							memcpy((char*)mxGetPr(M), (char*)&mNt, 1 * sizeof(double));
							engPutVariable(ep, "TX", T);
							engPutVariable(ep, "TY", T2);
							engPutVariable(ep, "TZ", T3);
							engPutVariable(ep, "NT", M);
							if (flag == false)
							{
								flag = true;
								engEvalString(ep, "h=plot3(TX,TY,TZ,'ko'); view(135,30); grid on;");
								az::mea::cmea::ePlot(ep, problem, pEA->Population().P().FSize());
							}
							else
							{
								engEvalString(ep, "set(h, 'XData',TX     ,'YData',TY  ,'ZData',TZ  )");
							}
						}
					}
					//engEvalString(ep,"hold off");
					//engClose(ep);
				}
                az::mea::cmea::ShowProgress(ir, generation, ic, max_gen);

			}
			// save data
			//std::stringstream spfs, spss;
			//spfs << "data/" << method.c_str() << "_" << problem.c_str() << "_" << ir << ".py";
			////std::cout<<spfs.str()<<std::endl;
			//f1.open(spfs.str().c_str());
			//f1 << std::scientific << std::setprecision(8);
			//for (i = 0; i<ipf; i++)
			//{
			//	for (j = 0; j<pEA->P().FSize(); j++) f1 << PFs[i][j] << "\t";
			//	f1 << std::endl;
			//}
			//f1.close();

			//spss << "data/" << method.c_str() << "_" << problem.c_str() << "_" << ir << ".px";
			//f2.open(spss.str().c_str());
			////std::cout<<ss0.str()<<std::endl;
			//f2 << std::scientific << std::setprecision(8);
			//for (i = 0; i<ipf; i++)
			//{
			//	for (j = 0; j<xdim; j++) f2 << PSs[i][j] << "\t";
			//	f2 << std::endl;
			//}
			//f2.close();

			//PFs.clear();
			//PSs.clear();

			// save final population
			std::stringstream spf;
			spf << "PF/" << method.c_str() << "_" << problem.c_str() << "_" << ir << ".dat";
			pf.open(spf.str().c_str());
			pf << std::setiosflags(std::ios::fixed) << std::setprecision(8);
			if (pEA->Population().FN() != 0) 
			for (i = 0; i<pEA->Population().Size(); i++)
			{
				for (j = 0; j<pEA->P().FSize(); j++) pf << pEA->Population()[i].F(j) << "\t";
				pf << std::endl;
			}
			pf.close();

			if (!DebugModel)
			{
				std::vector<vector<double>> ePf;
				for (unsigned int ei = 0; ei < pEA->Population().Size(); ei++) if (pEA->Population().In(ei).IsFeasible())
				{
					ePf.push_back(pEA->Population().In(ei).F());
				}
				if (pEA->Population().FN() != 0)	
				{
					pEval->FinalResults(ePf);
				}
			}
		}
		finish = clock();
		std::cout << "\r\t\t\t\t\t\t\t\t\r  ->> Completed !	Total Time: " << std::fixed << setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << std::endl;

		if (!DebugModel)	{ pEval->OutputTotal(); }
		delete pEval;
		delete pEA;
		pro_index++;
	}
	getchar();
	return 1;
}