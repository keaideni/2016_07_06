#include "Parameter.h"
#include "OP.h"
#include "SuperEnergy.h"
#include "DMRGP.h"
#include "physics.h"



int OP::Max;
int main()
{
	Parameter para;
	para.read();
	
	OP::Max = para.ParticleNo + 1;

	
	//para.D = 200;
	





	DMRGP DMRG(para);
	//DMRG.fwave.show();

	std::ofstream Fdata("./result/data", std::ios_base::out | std::ios_base::app);
	CalcuCorr(DMRG.OrbitalM, DMRG.fwave, Fdata);
	calcustructure(DMRG.fwave, DMRG.OrbitalM, para.ParticleNo, Fdata);
	Fdata.close();


	calcudensity(DMRG.OrbitalM, DMRG.fwave, para.ParticleNo);

	
}