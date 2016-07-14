#include "QWave.h"
#include "Corr.h"
#include "DMRGP.h"
//================function define===================================================
void CalcuCorr(const int& OrbitalM, const QWave& fwave, std::ofstream& Fdata);
double CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag, const QWave& fwave);
double CacuCorr(const OP& corrc, const OP& corrcdag, const QWave& fwave);


//=================calculate two types of correlation function========================
void CalcuCorr(const int& OrbitalM, const QWave& fwave, std::ofstream& Fdata)
{
        std::ofstream outfile1, outfile11;//to save the normalized correlation.
        std::ofstream outfile2, outfile22;//to save the not normalized correlation.
        if(OrbitalM % 2 == 1)
        {
                                        
                outfile1.open("./result/resonatorN.txt");//to store the type of particle besides the OrbitalM.
                outfile2.open("./result/resonator.txt");//to store the type of particle besides the OrbitalM.
                
                outfile11.open("./result/qubitN.txt");//to store the type of particle next besides the OrbitalM.
                outfile22.open("./result/qubit.txt");//to store the type of particles next besides the orbitalM.
        }else
        {
                                        
                outfile1.open("./result/qubitN.txt");//to store the type of particle besides the OrbitalM.
                outfile2.open("./result/qubit.txt");//to store the type of particles besides the orbitalM.
                
                outfile11.open("./result/resonatorN.txt");//to store the type of particle next besides the OrbitalM.
                outfile22.open("./result/resonator");//to store the type of particles next besides the orbitalM.
        }



        int i(OrbitalM + 1);//this label for Sigma;
        int j(OrbitalM - 1);//this label for Sigmadag and N;
        int fflag(1);

        double CorrLenth(0);
        double CorrSum(0);
        double correlation(0);

        Corr corr, corrdag, corrn;
        while(true)
        {
                                        
                                        

                corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                                        //corrn.show();Sys.SubSysEye.show();
                //corr.show();
                correlation = CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO, fwave);


                int distence(i - j);

                std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                outfile1<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;


                correlation = CacuCorr(corr.CorrO, corrdag.CorrO, fwave);
                std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                outfile2<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;

                CorrLenth += pow((i-j)/2, 2)*correlation;
                
                CorrSum +=correlation;
        

                if(fflag == 1)
                {
                        i += 2;
                }else
                {
                        j-= 2;
                }
                if((i>(OrbitalM-1) * 2)||(j<0))break;
                fflag *= -1;


        }
        
        CorrLenth = CorrLenth/CorrSum;
        if(OrbitalM %2 == 0)
        {
                Fdata << "the Qubit correlation = " <<CorrLenth<<std::endl;
                std::cout << "the Qubit correlated length = " <<CorrLenth<<std::endl;
        }else
        {
                Fdata << "the Resonator correlation = " <<CorrLenth<<std::endl;
                std::cout << "the Resonator correlated length = " <<CorrLenth<<std::endl;
        }



        i=(OrbitalM + 2);//this label for Sigma;
        j=(OrbitalM - 2);//this label for Sigmadag and N;
        fflag=1;
        CorrLenth = 0;
        CorrSum = 0;
        while(true)
        {
                                        
                                        

                corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                //corrn.show();Sys.SubSysEye.show();

                correlation = CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO, fwave);

                int distence(i - j);

                std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                outfile11<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;


                correlation = CacuCorr(corr.CorrO, corrdag.CorrO, fwave);
                std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                outfile22<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                CorrLenth += pow((i-j)/2, 2)*correlation;
                CorrSum += correlation;

                if(fflag == 1)
                {
                        i += 2;
                }else
                {
                        j-= 2;
                }
                if((i>(OrbitalM-1)*2)||(j<0))break;
                fflag *= -1;


        }

        CorrLenth /= CorrSum;

        if(OrbitalM %2 == 0)
        {
                Fdata << "the Resonator correlation = " <<CorrLenth<<std::endl;
                std::cout << "the Resonator correlated length = " <<CorrLenth<<std::endl;
        }else
        {
                Fdata << "the Qubit correlation = " <<CorrLenth<<std::endl;
                std::cout << "the Qubit correlated length = " <<CorrLenth<<std::endl;
        }

        outfile1.close();
        outfile2.close();
        outfile11.close();
        outfile22.close();
}
//=================================================================================================================


//=======================the first correlation which is normalized===========================================
double CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag, const QWave& fwave)
{
        

        QWave wave1, wave2;
        std::vector<double> f, f1;
        fwave.Wave2f(f);
        QWave ffwave(fwave);
        double number(0);
        double correlation(0);
        wave1.clear();
        wave1.OSWave2New(corrn, fwave);
        ffwave.initial(wave1);
        ffwave.Wave2f(f1);

        for(int i = 0; i < f.size(); ++i)
        {
                number += f[i]*f1[i];
        }

        wave1.clear();
        wave1.OSWave2New(corrcdag, fwave);
        wave2.clear();
        wave2.OEWave2New(corrc, wave1);
        ffwave.initial(wave2);
        ffwave.Wave2f(f1);
        correlation = 0;
        for(int i = 0; i < f.size(); ++i)
        {
                correlation += f[i]*f1[i];
        }
        correlation /= number;
        //std::cout<<correlation<<std::endl;
        
        return correlation;
        //std::cout<<number<<std::endl;


}
//=====================================================================================================

//================the second correlation which isn't normalized========================================
double CacuCorr(const OP& corrc, const OP& corrcdag, const QWave& fwave)
{
        QWave wave1, wave2;
        std::vector<double> f, f1;
        fwave.Wave2f(f);
        QWave ffwave(fwave);
        
        double correlation(0);

        wave1.clear();
        wave1.OSWave2New(corrcdag, fwave);
        wave2.clear();
        wave2.OEWave2New(corrc, wave1);
        ffwave.initial(wave2);
        ffwave.Wave2f(f1);
        correlation = 0;
        for(int i = 0; i < f.size(); ++i)
        {
                correlation += f[i]*f1[i];
        }
        
        
        return correlation;
        //std::cout<<number<<std::endl;


}
//====================used to calculate the correlation length which is defined directly================



//===============calculate the density========================================================
