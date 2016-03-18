#include <map>
#include <string>
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
using namespace std;
using namespace edm;
using namespace reco;

class MipCalibrator
{
public:
    MipCalibrator();
    
    float GetAct_fC2MIP(std::string type){return fC2MIP_act[type]; }
    float GetAct_MIP2keV(std::string type){return dEdX_act[type]; }
    float GetAbs_MIP2keV(int layer) {return dEdX_abs[layer]; }
    float GetAbs_weightkeV(int layer) {return 1.e03 * dEdX_absW[layer]; }

private:

    std::map<std::string, float> fC2MIP_act;
    std::map<std::string, float> dEdX_act;
    std::map<int, float> dEdX_abs;
    std::map<int, float> dEdX_absW;
};

MipCalibrator::MipCalibrator()
{
    fC2MIP_act["100"] = 1./1.43;   // fC2keV
    fC2MIP_act["200"] = 1./2.63;
    fC2MIP_act["300"] = 1./4.03;
    //    fC2MIP_act["BH"] = ; 

    dEdX_act["100"] = 27.55; //keV
    dEdX_act["200"] = 55.1; //keV
    dEdX_act["300"] = 85.0; //keV
    dEdX_act["BH"] = 633.7;

    dEdX_abs[1] = 7.298; //MeV
    dEdX_abs[2] = 9.908;
    dEdX_abs[3] = 6.227;
    dEdX_abs[4] = 9.908;
    dEdX_abs[5] = 6.227;
    dEdX_abs[6] = 9.908;
    dEdX_abs[7] = 6.227;
    dEdX_abs[8] = 9.908;
    dEdX_abs[9] = 6.227;
    dEdX_abs[10] = 9.908;
    dEdX_abs[11] = 7.995;
    dEdX_abs[12] = 12.275;
    dEdX_abs[13] = 7.995;
    dEdX_abs[14] = 12.275;
    dEdX_abs[15] = 7.995;
    dEdX_abs[16] = 12.275;
    dEdX_abs[17] = 7.995;
    dEdX_abs[18] = 12.275;
    dEdX_abs[19] = 7.995;
    dEdX_abs[20] = 12.275;
    dEdX_abs[21] = 11.089;
    dEdX_abs[22] = 16.219;
    dEdX_abs[23] = 11.089;
    dEdX_abs[24] = 16.219;
    dEdX_abs[25] = 11.089;
    dEdX_abs[26] = 16.219;
    dEdX_abs[27] = 11.089;
    dEdX_abs[28] = 16.219;
    dEdX_abs[29] = 60.182;
    dEdX_abs[30] = 49.871;
    dEdX_abs[31] = 49.871;
    dEdX_abs[32] = 49.871;
    dEdX_abs[33] = 49.871;
    dEdX_abs[34] = 49.871;
    dEdX_abs[35] = 49.871;
    dEdX_abs[36] = 49.871;
    dEdX_abs[37] = 49.871;
    dEdX_abs[38] = 49.871;
    dEdX_abs[39] = 49.871;
    dEdX_abs[40] = 49.871;
    dEdX_abs[41] = 74.139;
    dEdX_abs[42] = 92.196;
    dEdX_abs[43] = 92.196;
    dEdX_abs[44] = 92.196;
    dEdX_abs[45] = 92.196;
    dEdX_abs[46] = 92.196;
    dEdX_abs[47] = 92.196;
    dEdX_abs[48] = 92.196;
    dEdX_abs[49] = 92.196;
    dEdX_abs[50] = 92.196;
    dEdX_abs[51] = 92.196;
    dEdX_abs[52] = 92.196;
    // //convecrt in GeV
    // for(unsigned int iC=0; iC<dEdX_abs.size(); ++iC) dEdX_abs[iC] = dEdX_abs[iC]/1000.;


    dEdX_absW[1] = 8.603; // Mev
    dEdX_absW[2] = 8.0675; 
    dEdX_absW[3] = 8.0675; 
    dEdX_absW[4] = 8.0675; 
    dEdX_absW[5] = 8.0675; 
    dEdX_absW[6] = 8.0675; 
    dEdX_absW[7] = 8.0675; 
    dEdX_absW[8] = 8.0675; 
    dEdX_absW[9] = 8.0675; 
    dEdX_absW[10] = 8.9515; 
    dEdX_absW[11] = 10.135; 
    dEdX_absW[12] = 10.135; 
    dEdX_absW[13] = 10.135; 
    dEdX_absW[14] = 10.135; 
    dEdX_absW[15] = 10.135; 
    dEdX_absW[16] = 10.135; 
    dEdX_absW[17] = 10.135; 
    dEdX_absW[18] = 10.135; 
    dEdX_absW[19] = 10.135; 
    dEdX_absW[20] = 11.682; 
    dEdX_absW[21] = 13.654; 
    dEdX_absW[22] = 13.654; 
    dEdX_absW[23] = 13.654; 
    dEdX_absW[24] = 13.654; 
    dEdX_absW[25] = 13.654; 
    dEdX_absW[26] = 13.654; 
    dEdX_absW[27] = 13.654; 
    dEdX_absW[28] = 38.2005; 
    dEdX_absW[29] = 55.0265; 
    dEdX_absW[30] = 49.871; 
    dEdX_absW[31] = 49.871; 
    dEdX_absW[32] = 49.871; 
    dEdX_absW[33] = 49.871; 
    dEdX_absW[34] = 49.871; 
    dEdX_absW[35] = 49.871; 
    dEdX_absW[36] = 49.871; 
    dEdX_absW[37] = 49.871; 
    dEdX_absW[38] = 49.871; 
    dEdX_absW[39] = 49.871; 
    dEdX_absW[40] = 62.005; 
    dEdX_absW[41] = 83.1675; 
    dEdX_absW[42] = 92.196; 
    dEdX_absW[43] = 92.196; 
    dEdX_absW[44] = 92.196; 
    dEdX_absW[45] = 92.196; 
    dEdX_absW[46] = 92.196; 
    dEdX_absW[47] = 92.196; 
    dEdX_absW[48] = 92.196; 
    dEdX_absW[49] = 92.196; 
    dEdX_absW[50] = 92.196; 
    dEdX_absW[51] = 92.196; 
    dEdX_absW[52] = 46.098; 

}

//DEFINE_FWK_MODULE( MipCalibrator );

