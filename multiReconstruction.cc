#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TString.h"
#include "data_model.hh"
#include "RecHit.hh"
#include "Hit2Cluster.hh"
#include "Cluster.hh"
#include "Shower.hh"
#include "Neighbor.hh"
#include "Parameter.hh"
#include "SeedFinder.hh"
#include "Cluster2Shower.hh"
#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <map>
#include <vector>
#include "TChain.h"

using namespace std;

bool isTextFile(const string &fileName)
{
  size_t dotPos = fileName.find_last_of('.');
  if (dotPos == std::string::npos)
    return false;

  std::string extension = fileName.substr(dotPos + 1);
  if (extension == "txt")
    return true;
  else
    return false;
}

int main(int argc, char const *argv[])
{
  // 单模板pipeline对25个通道的时幅修正系数（使用数据修正）
  // double calibration_factor[25][25] = {{1.01278, 1.01317, 1.01314, 1.01188, 1.01224, 1.01202, 1.01194, 1.01133, 1.01088, 1.01159, 1.01082, 1.01079, 1.01084, 1.01083, 1.01151, 1.0116, 1.00985, 1.01183, 1.01264, 1.01293, 1.01115, 1.011, 1.01267, 1.01253, 1.01506}, {1.0106, 1.01108, 1.01063, 1.01032, 1.01007, 1.00985, 1.0104, 1.00941, 1.01018, 1.00965, 1.0097, 1.00948, 1.00934, 1.00898, 1.00934, 1.00959, 1.01007, 1.00996, 1.00997, 1.00961, 1.01089, 1.01074, 1.01135, 1.01215, 1.0126}, {1.00969, 1.00935, 1.00901, 1.00894, 1.00867, 1.00818, 1.00807, 1.00822, 1.00792, 1.00818, 1.00858, 1.00808, 1.00743, 1.00814, 1.00721, 1.00786, 1.00773, 1.00728, 1.00858, 1.00785, 1.00823, 1.00851, 1.00911, 1.00914, 1.0102}, {1.00622, 1.00493, 1.00514, 1.00487, 1.00376, 1.00368, 1.00399, 1.0043, 1.00296, 1.00241, 1.00314, 1.00377, 1.00421, 1.00335, 1.0035, 1.00389, 1.00298, 1.00347, 1.00284, 1.00382, 1.00367, 1.00393, 1.0049, 1.00472, 1.00465}, {1.00849, 1.00825, 1.00773, 1.00698, 1.00764, 1.00688, 1.00697, 1.0066, 1.00678, 1.00772, 1.00736, 1.00727, 1.00789, 1.00651, 1.00653, 1.00615, 1.00555, 1.00626, 1.00757, 1.00668, 1.00672, 1.00881, 1.0074, 1.00623, 1.00894}, {1.00872, 1.00883, 1.00728, 1.00807, 1.00809, 1.008, 1.00786, 1.00674, 1.00652, 1.00713, 1.00727, 1.00645, 1.00719, 1.00752, 1.00778, 1.00638, 1.00655, 1.00748, 1.00694, 1.00764, 1.00754, 1.0071, 1.00785, 1.00851, 1.00917}, {1.00926, 1.00924, 1.00875, 1.00884, 1.00848, 1.00821, 1.00819, 1.00822, 1.00786, 1.00773, 1.00764, 1.00749, 1.0076, 1.00766, 1.00749, 1.00752, 1.00749, 1.00785, 1.00796, 1.00801, 1.0084, 1.00892, 1.00922, 1.00949, 1.01018}, {1.00723, 1.00689, 1.00668, 1.00642, 1.00611, 1.00593, 1.006, 1.00578, 1.00574, 1.00527, 1.00526, 1.00524, 1.00513, 1.00512, 1.00501, 1.00493, 1.0051, 1.00516, 1.00537, 1.00567, 1.00554, 1.00615, 1.00642, 1.0068, 1.00729}, {0.994468, 0.994229, 0.993808, 0.993719, 0.993425, 0.993269, 0.993241, 0.992649, 0.992561, 0.992192, 0.992475, 0.992136, 0.992015, 0.991855, 0.992123, 0.991915, 0.992158, 0.992099, 0.992107, 0.992148, 0.991933, 0.992536, 0.992556, 0.992555, 0.992699}, {1.01084, 1.0103, 1.00935, 1.00919, 1.01066, 1.00985, 1.01037, 1.00993, 1.00844, 1.00952, 1.00863, 1.00925, 1.00882, 1.00873, 1.00934, 1.00828, 1.0096, 1.00968, 1.00925, 1.00997, 1.01046, 1.01065, 1.01013, 1.01132, 1.01159}, {1.01297, 1.01285, 1.01276, 1.01258, 1.01149, 1.01233, 1.01158, 1.01111, 1.01162, 1.01154, 1.01045, 1.01137, 1.0115, 1.01019, 1.01068, 1.01066, 1.01207, 1.01133, 1.01162, 1.01125, 1.01305, 1.01348, 1.01414, 1.01355, 1.01349}, {1.01743, 1.01721, 1.01703, 1.01679, 1.01676, 1.01654, 1.01633, 1.01606, 1.01623, 1.01585, 1.01599, 1.01556, 1.01566, 1.01546, 1.01578, 1.01591, 1.01608, 1.01615, 1.01657, 1.01676, 1.01725, 1.01771, 1.01821, 1.01866, 1.01897}, {1.00212, 1.00199, 1.00195, 1.00165, 1.0016, 1.00139, 1.00126, 1.00124, 1.00102, 1.00086, 1.00071, 1.00048, 1.00034, 1.00023, 1, 0.999869, 0.999806, 0.999768, 0.999731, 0.999755, 1.00001, 1.00022, 1.00066, 1.00114, 1.00166}, {1.00795, 1.00787, 1.00759, 1.00757, 1.00736, 1.00683, 1.00691, 1.00686, 1.00655, 1.00648, 1.00622, 1.0062, 1.00619, 1.00595, 1.006, 1.00625, 1.00601, 1.00593, 1.00635, 1.00673, 1.00675, 1.00724, 1.00756, 1.00805, 1.00851}, {1.01643, 1.01566, 1.01569, 1.0155, 1.01456, 1.01606, 1.01506, 1.01495, 1.01519, 1.01482, 1.01336, 1.01455, 1.01441, 1.01546, 1.01465, 1.01431, 1.01499, 1.01552, 1.01583, 1.01501, 1.01586, 1.01578, 1.01704, 1.01672, 1.01848}, {1.01273, 1.01241, 1.01141, 1.01128, 1.01121, 1.01181, 1.01205, 1.01131, 1.01013, 1.01081, 1.01037, 1.0092, 1.00928, 1.01084, 1.01046, 1.0095, 1.0107, 1.01099, 1.01095, 1.01234, 1.01205, 1.0115, 1.01209, 1.01243, 1.01308}, {1.01084, 1.01079, 1.01052, 1.01003, 1.00986, 1.00968, 1.00953, 1.0093, 1.00923, 1.00938, 1.00905, 1.00887, 1.00876, 1.00886, 1.00902, 1.00881, 1.00899, 1.00924, 1.00916, 1.00913, 1.01009, 1.01017, 1.01054, 1.01091, 1.01134}, {1.00531, 1.00511, 1.00486, 1.00477, 1.0044, 1.00416, 1.00409, 1.00396, 1.00387, 1.00343, 1.00331, 1.00341, 1.00311, 1.0032, 1.00319, 1.00313, 1.00306, 1.00335, 1.00362, 1.00379, 1.00377, 1.00415, 1.00447, 1.0047, 1.00518}, {1.0084, 1.00798, 1.00821, 1.00785, 1.00788, 1.00739, 1.00749, 1.00734, 1.00714, 1.00669, 1.00709, 1.00686, 1.00657, 1.00658, 1.00661, 1.00655, 1.0068, 1.00645, 1.00694, 1.00752, 1.0075, 1.00791, 1.00813, 1.00829, 1.00863}, {1.01306, 1.01291, 1.01183, 1.01208, 1.0133, 1.01237, 1.01289, 1.01164, 1.01164, 1.01171, 1.01302, 1.01258, 1.01245, 1.01156, 1.01226, 1.01215, 1.01057, 1.01235, 1.01231, 1.01173, 1.01288, 1.01215, 1.01376, 1.0148, 1.01421}, {1.01258, 1.01368, 1.01391, 1.01334, 1.01239, 1.01202, 1.01349, 1.01205, 1.01178, 1.01231, 1.01121, 1.01193, 1.01276, 1.01099, 1.01154, 1.01179, 1.01154, 1.01179, 1.01303, 1.01139, 1.01288, 1.01343, 1.01331, 1.01473, 1.01317}, {1.01591, 1.01623, 1.01638, 1.0159, 1.01451, 1.01459, 1.01487, 1.01462, 1.01418, 1.01445, 1.01459, 1.01566, 1.0144, 1.01421, 1.01419, 1.01491, 1.01463, 1.01409, 1.01617, 1.01568, 1.01607, 1.01645, 1.01702, 1.01613, 1.01702}, {1.00193, 1.00215, 1.00246, 1.00087, 1.00183, 1.00186, 1.00128, 1.00088, 1.00011, 1.00022, 1.0008, 0.999774, 0.999952, 1.00054, 0.999562, 1.00072, 1.00017, 1.00035, 1.00007, 1.00125, 1.00144, 1.00073, 1.00152, 1.00167, 1.0018}, {1.0104, 1.00935, 1.00993, 1.00853, 1.00919, 1.00851, 1.00838, 1.00792, 1.00742, 1.00856, 1.00745, 1.0081, 1.00856, 1.00641, 1.00755, 1.00806, 1.00734, 1.00802, 1.00835, 1.00873, 1.00892, 1.00848, 1.00938, 1.00825, 1.00968}, {1.01025, 1.01046, 1.00811, 1.00928, 1.01056, 1.00968, 1.00947, 1.00823, 1.00906, 1.00914, 1.00881, 1.00934, 1.00728, 1.00761, 1.00846, 1.00859, 1.00745, 1.00898, 1.00912, 1.00918, 1.0103, 1.00941, 1.00933, 1.00947, 1.01023}};
  // 理论曲线计算的修正系数
  double correct[1250] = {0.997872945, 0.997879794, 0.997886632, 0.997893459, 0.997900275, 0.997907079, 0.997913872, 0.997920654, 0.997927425, 0.997934185, 0.997940934, 0.997947671, 0.997954398, 0.997961113, 0.997967817, 0.99797451, 0.997981191, 0.997987862, 0.997994521, 0.998001169, 0.998007806, 0.998014432, 0.998021047, 0.998027651, 0.998034243, 0.998040824, 0.998047395, 0.998053954, 0.998060502, 0.998067038, 0.998073564, 0.998080078, 0.998086581, 0.998093074, 0.998099555, 0.998106024, 0.998112483, 0.998118931, 0.998125367, 0.998131792, 0.998138206, 0.998144609, 0.998151001, 0.998157382, 0.998163751, 0.99817011, 0.998176457, 0.998182793, 0.998189118, 0.998195432, 0.998201735, 0.998208026, 0.998214307, 0.998220576, 0.998226834, 0.998233081, 0.998239317, 0.998245542, 0.998251755, 0.998257958, 0.998264149, 0.99827033, 0.998276499, 0.998282657, 0.998288803, 0.998294939, 0.998301064, 0.998307177, 0.998313279, 0.998319371, 0.998325451, 0.99833152, 0.998337577, 0.998343624, 0.998349659, 0.998355684, 0.998361697, 0.998367699, 0.99837369, 0.99837967, 0.998385639, 0.998391597, 0.998397543, 0.998403479, 0.998409403, 0.998415316, 0.998421218, 0.998427109, 0.998432989, 0.998438858, 0.998444715, 0.998450562, 0.998456397, 0.998462221, 0.998468035, 0.998473837, 0.998479627, 0.998485407, 0.998491176, 0.998496933, 0.99850268, 0.998508415, 0.998514139, 0.998519852, 0.998525554, 0.998531245, 0.998536925, 0.998542594, 0.998548251, 0.998553898, 0.998559533, 0.998565157, 0.998570771, 0.998576373, 0.998581963, 0.998587543, 0.998593112, 0.99859867, 0.998604216, 0.998609752, 0.998615276, 0.998620789, 0.998626291, 0.998631782, 0.998637262, 0.998642731, 0.998648188, 0.998653635, 0.998659071, 0.998664495, 0.998669908, 0.998675311, 0.998680702, 0.998686082, 0.998691451, 0.998696808, 0.998702155, 0.998707491, 0.998712815, 0.998718129, 0.998723431, 0.998728723, 0.998734003, 0.998739272, 0.99874453, 0.998749777, 0.998755013, 0.998760237, 0.998765451, 0.998770654, 0.998775845, 0.998781026, 0.998786195, 0.998791353, 0.9987965, 0.998801636, 0.998806761, 0.998811875, 0.998816978, 0.99882207, 0.998827151, 0.99883222, 0.998837279, 0.998842326, 0.998847363, 0.998852388, 0.998857402, 0.998862405, 0.998867397, 0.998872378, 0.998877348, 0.998882307, 0.998887255, 0.998892192, 0.998897117, 0.998902032, 0.998906935, 0.998911828, 0.998916709, 0.998921579, 0.998926438, 0.998931287, 0.998936124, 0.99894095, 0.998945765, 0.998950569, 0.998955361, 0.998960143, 0.998964914, 0.998969673, 0.998974422, 0.99897916, 0.998983886, 0.998988601, 0.998993306, 0.998997999, 0.999002681, 0.999007352, 0.999012012, 0.999016661, 0.999021299, 0.999025926, 0.999030542, 0.999035147, 0.999039741, 0.999044323, 0.999048895, 0.999053456, 0.999058005, 0.999062544, 0.999067071, 0.999071588, 0.999076093, 0.999080587, 0.99908507, 0.999089543, 0.999094004, 0.999098454, 0.999102893, 0.999107321, 0.999111738, 0.999116144, 0.999120539, 0.999124923, 0.999129296, 0.999133657, 0.999138008, 0.999142348, 0.999146676, 0.999150994, 0.999155301, 0.999159596, 0.999163881, 0.999168154, 0.999172416, 0.999176668, 0.999180908, 0.999185137, 0.999189356, 0.999193563, 0.999197759, 0.999201944, 0.999206119, 0.999210282, 0.999214434, 0.999218575, 0.999222705, 0.999226824, 0.999230932, 0.999235029, 0.999239115, 0.99924319, 0.999247254, 0.999251307, 0.999255348, 0.999259379, 0.999263399, 0.999267408, 0.999271405, 0.999275392, 0.999279368, 0.999283333, 0.999287286, 0.999291229, 0.999295161, 0.999299081, 0.999302991, 0.999306889, 0.999310777, 0.999314654, 0.999318519, 0.999322374, 0.999326217, 0.99933005, 0.999333871, 0.999337682, 0.999341481, 0.99934527, 0.999349047, 0.999352814, 0.999356569, 0.999360313, 0.999364047, 0.999367769, 0.999371481, 0.999375181, 0.999378871, 0.999382549, 0.999386217, 0.999389873, 0.999393518, 0.999397153, 0.999400776, 0.999404389, 0.99940799, 0.99941158, 0.99941516, 0.999418728, 0.999422286, 0.999425832, 0.999429368, 0.999432892, 0.999436406, 0.999439908, 0.9994434, 0.99944688, 0.99945035, 0.999453808, 0.999457256, 0.999460692, 0.999464118, 0.999467532, 0.999470936, 0.999474328, 0.99947771, 0.999481081, 0.99948444, 0.999487789, 0.999491127, 0.999494453, 0.999497769, 0.999501074, 0.999504368, 0.99950765, 0.999510922, 0.999514183, 0.999517433, 0.999520672, 0.999523899, 0.999527116, 0.999530322, 0.999533517, 0.999536701, 0.999539874, 0.999543036, 0.999546187, 0.999549328, 0.999552457, 0.999555575, 0.999558682, 0.999561778, 0.999564864, 0.999567938, 0.999571001, 0.999574054, 0.999577095, 0.999580125, 0.999583145, 0.999586153, 0.999589151, 0.999592138, 0.999595113, 0.999598078, 0.999601032, 0.999603974, 0.999606906, 0.999609827, 0.999612737, 0.999615636, 0.999618524, 0.999621401, 0.999624267, 0.999627122, 0.999629966, 0.999632799, 0.999635622, 0.999638433, 0.999641233, 0.999644023, 0.999646801, 0.999649569, 0.999652325, 0.999655071, 0.999657806, 0.99966053, 0.999663242, 0.999665944, 0.999668635, 0.999671315, 0.999673984, 0.999676642, 0.999679289, 0.999681926, 0.999684551, 0.999687165, 0.999689769, 0.999692361, 0.999694943, 0.999697513, 0.999700073, 0.999702622, 0.999705159, 0.999707686, 0.999710202, 0.999712707, 0.999715201, 0.999717684, 0.999720157, 0.999722618, 0.999725068, 0.999727508, 0.999729936, 0.999732354, 0.999734761, 0.999737156, 0.999739541, 0.999741915, 0.999744278, 0.99974663, 0.999748971, 0.999751301, 0.999753621, 0.999755929, 0.999758227, 0.999760513, 0.999762789, 0.999765053, 0.999767307, 0.99976955, 0.999771782, 0.999774003, 0.999776213, 0.999778413, 0.999780601, 0.999782778, 0.999784945, 0.9997871, 0.999789245, 0.999791379, 0.999793502, 0.999795614, 0.999797715, 0.999799805, 0.999801884, 0.999803953, 0.99980601, 0.999808057, 0.999810093, 0.999812117, 0.999814131, 0.999816134, 0.999818126, 0.999820107, 0.999822078, 0.999824037, 0.999825986, 0.999827923, 0.99982985, 0.999831766, 0.999833671, 0.999835565, 0.999837448, 0.99983932, 0.999841182, 0.999843032, 0.999844872, 0.9998467, 0.999848518, 0.999850325, 0.999852121, 0.999853906, 0.999855681, 0.999857444, 0.999859197, 0.999860938, 0.999862669, 0.999864389, 0.999866098, 0.999867796, 0.999869484, 0.99987116, 0.999872825, 0.99987448, 0.999876124, 0.999877757, 0.999879379, 0.99988099, 0.99988259, 0.99988418, 0.999885758, 0.999887326, 0.999888883, 0.999890429, 0.999891964, 0.999893488, 0.999895001, 0.999896504, 0.999897995, 0.999899476, 0.999900946, 0.999902405, 0.999903853, 0.999905291, 0.999906717, 0.999908133, 0.999909538, 0.999910931, 0.999912314, 0.999913687, 0.999915048, 0.999916398, 0.999917738, 0.999919067, 0.999920385, 0.999921692, 0.999922988, 0.999924274, 0.999925548, 0.999926812, 0.999928065, 0.999929307, 0.999930538, 0.999931758, 0.999932967, 0.999934166, 0.999935354, 0.999936531, 0.999937697, 0.999938852, 0.999939997, 0.99994113, 0.999942253, 0.999943365, 0.999944466, 0.999945556, 0.999946636, 0.999947704, 0.999948762, 0.999949809, 0.999950845, 0.99995187, 0.999952885, 0.999953888, 0.999954881, 0.999955863, 0.999956834, 0.999957794, 0.999958744, 0.999959682, 0.99996061, 0.999961527, 0.999962433, 0.999963329, 0.999964213, 0.999965087, 0.99996595, 0.999966802, 0.999967643, 0.999968474, 0.999969293, 0.999970102, 0.9999709, 0.999971687, 0.999972464, 0.999973229, 0.999973984, 0.999974728, 0.999975461, 0.999976183, 0.999976895, 0.999977596, 0.999978285, 0.999978965, 0.999979633, 0.99998029, 0.999980937, 0.999981573, 0.999982198, 0.999982812, 0.999983416, 0.999984008, 0.99998459, 0.999985161, 0.999985722, 0.999986271, 0.99998681, 0.999987338, 0.999987855, 0.999988361, 0.999988857, 0.999989342, 0.999989816, 0.999990279, 0.999990731, 0.999991173, 0.999991604, 0.999992024, 0.999992433, 0.999992831, 0.999993219, 0.999993596, 0.999993962, 0.999994317, 0.999994662, 0.999994996, 0.999995319, 0.999995631, 0.999995932, 0.999996223, 0.999996503, 0.999996772, 0.99999703, 0.999997278, 0.999997514, 0.99999774, 0.999997956, 0.99999816, 0.999998354, 0.999998537, 0.999998709, 0.99999887, 0.999999021, 0.999999161, 0.99999929, 0.999999408, 0.999999516, 0.999999613, 0.999999699, 0.999999774, 0.999999839, 0.999999892, 0.999999935, 0.999999968, 0.999999989, 1, 1, 0.999999989, 0.999999968, 0.999999935, 0.999999892, 0.999999839, 0.999999774, 0.999999699, 0.999999613, 0.999999516, 0.999999409, 0.99999929, 0.999999161, 0.999999022, 0.999998871, 0.99999871, 0.999998538, 0.999998355, 0.999998162, 0.999997958, 0.999997743, 0.999997517, 0.999997281, 0.999997034, 0.999996776, 0.999996507, 0.999996228, 0.999995938, 0.999995637, 0.999995325, 0.999995003, 0.99999467, 0.999994327, 0.999993972, 0.999993607, 0.999993231, 0.999992844, 0.999992447, 0.999992039, 0.99999162, 0.999991191, 0.999990751, 0.9999903, 0.999989838, 0.999989366, 0.999988883, 0.999988389, 0.999987884, 0.999987369, 0.999986843, 0.999986306, 0.999985759, 0.999985201, 0.999984632, 0.999984053, 0.999983462, 0.999982861, 0.99998225, 0.999981628, 0.999980995, 0.999980351, 0.999979696, 0.999979031, 0.999978355, 0.999977669, 0.999976972, 0.999976264, 0.999975545, 0.999974816, 0.999974076, 0.999973325, 0.999972563, 0.999971791, 0.999971008, 0.999970215, 0.999969411, 0.999968596, 0.99996777, 0.999966934, 0.999966087, 0.999965229, 0.999964361, 0.999963482, 0.999962592, 0.999961692, 0.999960781, 0.999959859, 0.999958927, 0.999957984, 0.99995703, 0.999956065, 0.99995509, 0.999954104, 0.999953108, 0.999952101, 0.999951083, 0.999950054, 0.999949015, 0.999947965, 0.999946905, 0.999945833, 0.999944752, 0.999943659, 0.999942556, 0.999941442, 0.999940317, 0.999939182, 0.999938036, 0.99993688, 0.999935713, 0.999934535, 0.999933346, 0.999932147, 0.999930937, 0.999929717, 0.999928486, 0.999927244, 0.999925991, 0.999924728, 0.999923454, 0.99992217, 0.999920875, 0.999919569, 0.999918253, 0.999916926, 0.999915588, 0.99991424, 0.999912881, 0.999911511, 0.999910131, 0.99990874, 0.999907338, 0.999905926, 0.999904503, 0.99990307, 0.999901626, 0.999900171, 0.999898706, 0.99989723, 0.999895743, 0.999894246, 0.999892738, 0.999891219, 0.99988969, 0.99988815, 0.9998866, 0.999885039, 0.999883467, 0.999881885, 0.999880292, 0.999878688, 0.999877074, 0.999875449, 0.999873814, 0.999872168, 0.999870511, 0.999868844, 0.999867166, 0.999865477, 0.999863778, 0.999862068, 0.999860348, 0.999858617, 0.999856875, 0.999855123, 0.99985336, 0.999851586, 0.999849802, 0.999848008, 0.999846202, 0.999844386, 0.99984256, 0.999840723, 0.999838875, 0.999837017, 0.999835148, 0.999833268, 0.999831378, 0.999829477, 0.999827566, 0.999825644, 0.999823711, 0.999821768, 0.999819815, 0.99981785, 0.999815875, 0.99981389, 0.999811894, 0.999809887, 0.99980787, 0.999805842, 0.999803803, 0.999801754, 0.999799695, 0.999797624, 0.999795544, 0.999793452, 0.99979135, 0.999789238, 0.999787114, 0.999784981, 0.999782836, 0.999780681, 0.999778516, 0.99977634, 0.999774153, 0.999771956, 0.999769748, 0.99976753, 0.999765301, 0.999763062, 0.999760812, 0.999758551, 0.99975628, 0.999753998, 0.999751706, 0.999749403, 0.999747089, 0.999744765, 0.999742431, 0.999740085, 0.99973773, 0.999735363, 0.999732987, 0.999730599, 0.999728201, 0.999725793, 0.999723374, 0.999720944, 0.999718504, 0.999716053, 0.999713592, 0.99971112, 0.999708638, 0.999706145, 0.999703641, 0.999701127, 0.999698603, 0.999696067, 0.999693522, 0.999690966, 0.999688399, 0.999685822, 0.999683234, 0.999680635, 0.999678026, 0.999675407, 0.999672777, 0.999670136, 0.999667485, 0.999664824, 0.999662152, 0.999659469, 0.999656776, 0.999654072, 0.999651358, 0.999648633, 0.999645898, 0.999643152, 0.999640396, 0.999637629, 0.999634851, 0.999632063, 0.999629265, 0.999626456, 0.999623636, 0.999620806, 0.999617966, 0.999615115, 0.999612253, 0.999609381, 0.999606499, 0.999603605, 0.999600702, 0.999597788, 0.999594863, 0.999591928, 0.999588982, 0.999586026, 0.999583059, 0.999580082, 0.999577094, 0.999574096, 0.999571088, 0.999568068, 0.999565039, 0.999561998, 0.999558948, 0.999555887, 0.999552815, 0.999549733, 0.99954664, 0.999543537, 0.999540423, 0.999537299, 0.999534164, 0.999531019, 0.999527863, 0.999524697, 0.999521521, 0.999518334, 0.999515136, 0.999511928, 0.999508709, 0.99950548, 0.999502241, 0.999498991, 0.99949573, 0.999492459, 0.999489178, 0.999485886, 0.999482584, 0.999479271, 0.999475947, 0.999472613, 0.999469269, 0.999465914, 0.999462549, 0.999459173, 0.999455787, 0.999452391, 0.999448983, 0.999445566, 0.999442138, 0.999438699, 0.99943525, 0.999431791, 0.999428321, 0.999424841, 0.99942135, 0.999417849, 0.999414337, 0.999410815, 0.999407282, 0.999403739, 0.999400185, 0.999396621, 0.999393047, 0.999389462, 0.999385867, 0.999382261, 0.999378645, 0.999375018, 0.999371381, 0.999367733, 0.999364075, 0.999360407, 0.999356728, 0.999353039, 0.999349339, 0.999345629, 0.999341908, 0.999338177, 0.999334435, 0.999330683, 0.999326921, 0.999323148, 0.999319365, 0.999315571, 0.999311767, 0.999307952, 0.999304127, 0.999300292, 0.999296446, 0.99929259, 0.999288723, 0.999284846, 0.999280959, 0.999277061, 0.999273152, 0.999269233, 0.999265304, 0.999261365, 0.999257415, 0.999253454, 0.999249483, 0.999245502, 0.99924151, 0.999237508, 0.999233496, 0.999229473, 0.999225439, 0.999221396, 0.999217342, 0.999213277, 0.999209202, 0.999205117, 0.999201021, 0.999196915, 0.999192798, 0.999188671, 0.999184534, 0.999180386, 0.999176228, 0.99917206, 0.999167881, 0.999163691, 0.999159492, 0.999155281, 0.999151061, 0.99914683, 0.999142589, 0.999138337, 0.999134075, 0.999129803, 0.99912552, 0.999121227, 0.999116923, 0.999112609, 0.999108285, 0.99910395, 0.999099605, 0.999095249, 0.999090884, 0.999086507, 0.999082121, 0.999077724, 0.999073316, 0.999068899, 0.999064471, 0.999060032, 0.999055583, 0.999051124, 0.999046655, 0.999042175, 0.999037684, 0.999033184, 0.999028673, 0.999024151, 0.99901962, 0.999015078, 0.999010525, 0.999005962, 0.999001389, 0.998996806, 0.998992212, 0.998987608, 0.998982993, 0.998978368, 0.998973733, 0.998969087, 0.998964431, 0.998959765, 0.998955088, 0.998950401, 0.998945704, 0.998940996, 0.998936278, 0.99893155, 0.998926811, 0.998922062, 0.998917303, 0.998912533, 0.998907753, 0.998902963, 0.998898162, 0.998893351, 0.99888853, 0.998883698, 0.998878856, 0.998874003, 0.998869141, 0.998864268, 0.998859384, 0.998854491, 0.998849587, 0.998844672, 0.998839748, 0.998834813, 0.998829868, 0.998824912, 0.998819946, 0.99881497, 0.998809983, 0.998804986, 0.998799979, 0.998794962, 0.998789934, 0.998784896, 0.998779848, 0.998774789, 0.99876972, 0.99876464, 0.998759551, 0.998754451, 0.998749341, 0.99874422, 0.998739089, 0.998733948, 0.998728797, 0.998723635, 0.998718463, 0.998713281, 0.998708088, 0.998702885, 0.998697672, 0.998692448, 0.998687215, 0.998681971, 0.998676716, 0.998671452, 0.998666177, 0.998660892, 0.998655596, 0.99865029, 0.998644974, 0.998639648, 0.998634311, 0.998628964, 0.998623607, 0.99861824, 0.998612862, 0.998607474, 0.998602076, 0.998596667, 0.998591248, 0.998585819, 0.99858038, 0.99857493, 0.99856947, 0.998564, 0.99855852, 0.998553029, 0.998547528, 0.998542017, 0.998536496, 0.998530964, 0.998525422, 0.99851987, 0.998514307, 0.998508734, 0.998503151, 0.998497558, 0.998491955, 0.998486341, 0.998480717, 0.998475083, 0.998469438, 0.998463783, 0.998458118, 0.998452443, 0.998446758, 0.998441062, 0.998435356, 0.99842964, 0.998423913, 0.998418177, 0.99841243, 0.998406672, 0.998400905, 0.998395127, 0.99838934, 0.998383541, 0.998377733, 0.998371915, 0.998366086, 0.998360247, 0.998354398, 0.998348538, 0.998342668, 0.998336788, 0.998330898, 0.998324998, 0.998319087, 0.998313167, 0.998307236, 0.998301294, 0.998295343, 0.998289381, 0.998283409, 0.998277427, 0.998271435, 0.998265432, 0.99825942, 0.998253397, 0.998247364, 0.99824132, 0.998235267, 0.998229203, 0.998223129, 0.998217045, 0.998210951, 0.998204846, 0.998198731, 0.998192606, 0.998186471, 0.998180326, 0.99817417, 0.998168005, 0.998161829, 0.998155642, 0.998149446, 0.99814324, 0.998137023, 0.998130796, 0.998124559, 0.998118312, 0.998112054, 0.998105787, 0.998099509, 0.998093221, 0.998086923, 0.998080615, 0.998074296, 0.998067967, 0.998061628, 0.998055279, 0.99804892, 0.998042551, 0.998036171, 0.998029782, 0.998023382, 0.998016972, 0.998010551, 0.998004121, 0.997997681, 0.99799123, 0.997984769, 0.997978298, 0.997971817, 0.997965325, 0.997958824, 0.997952312, 0.99794579, 0.997939258};
  vector<string> datafiles;

  if (isTextFile(argv[1]))
  {
    ifstream filestream(argv[1]);
    if (filestream.is_open())
    {
      string filename;
      while (getline(filestream, filename))
        datafiles.push_back(filename);
    }
    else
      cout << "error in reading filelist  " << endl;
  }
  else
    datafiles.push_back(argv[1]);

  TChain *tree = new TChain("decode_data");

  for (size_t i = 0; i < datafiles.size(); i++)
    tree->AddFile(datafiles.at(i).data());

  int TriggerID;
  tree->SetBranchAddress("TriggerID", &TriggerID);
  vector<decode_data_col *> Hit(25);
  double TimeStamp[25];
  vector<vector<double> *> CoarseTime;
  vector<vector<double> *> FineTime;
  vector<vector<double> *> Amplitude;
  for (int i = 0; i < 25; i++)
  {
    CoarseTime.push_back(new vector<double>);
    FineTime.push_back(new vector<double>);
    Amplitude.push_back(new vector<double>);
  }
  string Name[25];
  for (int i = 0; i < 25; i++)
  {
    Name[i] = to_string(i / 5 + 1) + "_" + to_string(i % 5 + 1);
  }
  for (int i = 0; i < 25; i++)
    Hit[i] = new decode_data_col();
  for (int i = 0; i < 25; i++)
  {
    string name = "Hit_" + Name[i];
    string timeName = name + "_TimeStamp";
    string coarseName = name + "_CoarseTime";
    string fineName = name + "_FineTime";
    string ampName = name + "_Amplitude";
    tree->SetBranchAddress(name.data(), Hit[i]);
    tree->SetBranchAddress(timeName.data(), &TimeStamp[i]);
    tree->SetBranchAddress(coarseName.data(), &CoarseTime[i]);
    tree->SetBranchAddress(fineName.data(), &FineTime[i]);
    tree->SetBranchAddress(ampName.data(), &Amplitude[i]);
  }

  TString outputfile;
  if (argc >= 3)
    outputfile = argv[2];
  else
  {
    outputfile = "Reconstruction.root";
    cout << "Auto save file as Reconstruction.root..." << endl;
  }

  TFile *f = new TFile(outputfile.Data(), "recreate");

  int recmode = 1;
  if (argc == 4)
    recmode = atoi(argv[3]);

  int triggerID;
  vector<int> HitID;
  vector<double> HitEnergy;
  vector<double> HitTime;
  vector<int> ClusterID;
  vector<int> ClusterSize;
  vector<double> ClusterEnergy;
  vector<double> ClusterSMoment;
  vector<double> ClusterPosX;
  vector<double> ClusterPosY;
  vector<double> ClusterPosZ;
  vector<int> ShowerID;
  vector<double> ShowerE3x3;
  vector<double> ShowerE5x5;
  vector<double> ShowerEAll;
  vector<double> ShowerPosX3x3;
  vector<double> ShowerPosY3x3;
  vector<double> ShowerPosZ3x3;
  vector<double> ShowerPosX5x5;
  vector<double> ShowerPosY5x5;
  vector<double> ShowerPosZ5x5;
  vector<double> ShowerSMoment;
  vector<double> ShowerLatMoment;
  vector<double> ShowerA20Moment;
  vector<double> ShowerA42Moment;

  TTree *rec_data = new TTree("rec_data", "rec_data");
  rec_data->Branch("EventID", &triggerID, "triggerID/I");
  rec_data->Branch("HitID", &HitID);
  rec_data->Branch("HitEnergy", &HitEnergy);
  rec_data->Branch("HitTime", &HitTime);
  rec_data->Branch("ClusterID", &ClusterID);
  rec_data->Branch("ClusterSize", &ClusterSize);
  rec_data->Branch("ClusterEnergy", &ClusterEnergy);
  rec_data->Branch("ClusterSMoment", &ClusterSMoment);
  rec_data->Branch("ClusterPosX", &ClusterPosX);
  rec_data->Branch("ClusterPosY", &ClusterPosY);
  rec_data->Branch("ClusterPosZ", &ClusterPosZ);
  rec_data->Branch("ShowerID", &ShowerID);
  rec_data->Branch("ShowerE3x3", &ShowerE3x3);
  rec_data->Branch("ShowerE5x5", &ShowerE5x5);
  rec_data->Branch("ShowerEAll", &ShowerEAll);
  rec_data->Branch("ShowerPosX3x3", &ShowerPosX3x3);
  rec_data->Branch("ShowerPosY3x3", &ShowerPosY3x3);
  rec_data->Branch("ShowerPosZ3x3", &ShowerPosZ3x3);
  rec_data->Branch("ShowerPosX5x5", &ShowerPosX5x5);
  rec_data->Branch("ShowerPosY5x5", &ShowerPosY5x5);
  rec_data->Branch("ShowerPosZ5x5", &ShowerPosZ5x5);
  rec_data->Branch("ShowerSMoment", &ShowerSMoment);
  rec_data->Branch("ShowerLatMoment", &ShowerLatMoment);
  rec_data->Branch("ShowerA20Moment", &ShowerA20Moment);
  rec_data->Branch("ShowerA42Moment", &ShowerA42Moment);

  Parameter &Para = Parameter::GetInstance();

  TRandom3 r1;
  r1.SetSeed(time(0));

  const int npoint = 12;
  double hitenergy[npoint] = {0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.02, 0.024, 0.032, 0.085, 0.09, 0.1};
  double delta_time[npoint] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 10, 10};
  TGraph *timecut = new TGraph(npoint, hitenergy, delta_time);

  int nEntries = tree->GetEntries();
  int interval = nEntries / 20;
  for (int k = 0; k < nEntries; k++)
  {
    int progress = static_cast<float>(k + 1) / nEntries * 100;
    if ((k + 1) % interval == 0)
    {
      cout << "Progress: " << progress + 1 << "%\r" << endl;
      std::cout.flush();
    }

    tree->GetEntry(k);

    triggerID = TriggerID;

    HitID.clear();
    HitEnergy.clear();
    HitTime.clear();
    ClusterID.clear();
    ClusterSize.clear();
    ClusterEnergy.clear();
    ClusterSMoment.clear();
    ClusterPosX.clear();
    ClusterPosY.clear();
    ClusterPosZ.clear();
    ShowerID.clear();
    ShowerE3x3.clear();
    ShowerE5x5.clear();
    ShowerEAll.clear();
    ShowerPosX3x3.clear();
    ShowerPosY3x3.clear();
    ShowerPosZ3x3.clear();
    ShowerPosX5x5.clear();
    ShowerPosY5x5.clear();
    ShowerPosZ5x5.clear();
    ShowerSMoment.clear();
    ShowerLatMoment.clear();
    ShowerA20Moment.clear();
    ShowerA42Moment.clear();

    if (recmode == 0)
    {
      vector<int> SeedIDVec;
      map<int, RecHit> HitMap;
      map<int, Cluster> ClusterMap;
      map<int, Shower> ShowerMap;

      map<int, RecHit>::iterator Hiter;
      map<int, Cluster>::iterator Citer;
      map<int, Shower>::iterator Siter;

      SeedIDVec.clear();
      HitMap.clear();
      ClusterMap.clear();
      ShowerMap.clear();

      for (int i = 0; i < 25; i++)
      {
        int CrystalID = Hit[i]->CrystalID;
        double timestamp = TimeStamp[i];
        int index = -1;
        double maxAmp = -1;
        for (int j = 0; j < Amplitude[i]->size(); j++)
        {
          double coarseTime = CoarseTime[i]->at(j);
          double fineTime = FineTime[i]->at(j);
          double amp = Amplitude[i]->at(j);
          double time = coarseTime - timestamp;
          if (amp > maxAmp)
          {
            maxAmp = amp;
            index = j;
          }
        }
        if (index != -1)
        {
          RecHit rechit;
          rechit.Clear();
          rechit.setCrystalID(CrystalID);
          rechit.setFrontCenter(TVector3(Para.HitPosX(CrystalID), Para.HitPosY(CrystalID), Para.HitPosZ(CrystalID)));
          double coarseTime = CoarseTime[i]->at(index);
          double fineTime = FineTime[i]->at(index);
          double amp = Amplitude[i]->at(index);
          double time = coarseTime - timestamp;
          rechit.setEnergy(amp, Para.HGNoise(i), Para.HGPedestal(i), Para.HGMipPeak(i));
          rechit.setTime(time);
          if (rechit.Energy() > 0)
          {
            HitMap[CrystalID] = rechit;
            HitID.push_back(rechit.CrystalID());
            HitTime.push_back(rechit.Time());
            HitEnergy.push_back(rechit.Energy());
          }
        }
      }

      Hit2Cluster m_Hit2Cluster;
      m_Hit2Cluster.Convert(HitMap, ClusterMap);
      Cluster2Shower m_Cluster2Shower;
      m_Cluster2Shower.Convert(ClusterMap, ShowerMap);

      /*SeedFinder m_SeedFinder;
        Citer=ClusterMap.begin();
        m_SeedFinder.Seed(Citer->second,SeedIDVec);
        for(int i=0;i<SeedIDVec.size();i++) cout<<"Seed:"<<SeedIDVec.at(i)<<endl;*/

      double epsilon = 1e-10;
      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        // Calculate cluster position
        TVector3 possum;
        double etot = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          etot += it->second.Energy();
          TVector3 pos(it->second.FrontCenter());
          possum += pos * it->second.Energy();
        }
        if (etot > 0)
          possum = possum * (1 / etot);
        if (std::abs(possum.x() - (std::floor(possum.x() * 100000) / 100000)) < epsilon)
          Citer->second.setPosition(TVector3(std::floor(possum.x() * 100000) / 100000, std::floor(possum.y() * 100000) / 100000, std::floor(possum.z() * 100000) / 100000));
        else
          Citer->second.setPosition(possum);

        // Calculate cluster second moment
        double sum = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          TVector3 pos(it->second.FrontCenter());
          if (std::abs((pos - possum).Mag2()) < epsilon)
            sum += it->second.Energy() * 0;

          else
            sum += it->second.Energy() * ((pos - possum).Mag2());
        }
        if (etot > 0)
          sum /= etot;
        Citer->second.setSecondMoment(sum);
      }

      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        ClusterID.push_back(Citer->second.GetClusterID());
        ClusterSize.push_back(Citer->second.GetClusterSize());
        ClusterEnergy.push_back(Citer->second.GetClusterEnergy());
        ClusterSMoment.push_back(Citer->second.GetSecondMoment());
        ClusterPosX.push_back(Citer->second.GetPosition().x());
        ClusterPosY.push_back(Citer->second.GetPosition().y());
        ClusterPosZ.push_back(Citer->second.GetPosition().z());
      }

      for (Siter = ShowerMap.begin(); Siter != ShowerMap.end(); Siter++)
      {
        ShowerID.push_back(Siter->second.SeedID());
        ShowerE3x3.push_back(Siter->second.E3x3());
        ShowerE5x5.push_back(Siter->second.E5x5());
        ShowerEAll.push_back(Siter->second.EAll());
        ShowerPosX3x3.push_back(Siter->second.x3x3());
        ShowerPosY3x3.push_back(Siter->second.y3x3());
        ShowerPosZ3x3.push_back(Siter->second.z3x3());
        ShowerPosX5x5.push_back(Siter->second.x5x5());
        ShowerPosY5x5.push_back(Siter->second.y5x5());
        ShowerPosZ5x5.push_back(Siter->second.z5x5());
        ShowerSMoment.push_back(Siter->second.SecondMoment());
        ShowerLatMoment.push_back(Siter->second.LateralMoment());
        ShowerA20Moment.push_back(Siter->second.A20Moment());
        ShowerA42Moment.push_back(Siter->second.A42Moment());
      }
    }

    if (recmode == 1)
    {
      vector<int> SeedIDVec;
      map<int, RecHit> HitMap;
      map<int, Cluster> ClusterMap;
      map<int, Shower> ShowerMap;

      map<int, RecHit>::iterator Hiter;
      map<int, Cluster>::iterator Citer;
      map<int, Shower>::iterator Siter;

      SeedIDVec.clear();
      HitMap.clear();
      ClusterMap.clear();
      ShowerMap.clear();

      for (int i = 0; i < 25; i++)
      {
        int CrystalID = Hit[i]->CrystalID;
        double timestamp = TimeStamp[i];
        int index = -1;
        double maxAmp = -1;
                for (int j = 0; j < Amplitude[i]->size(); j++)
        {
          double coarseTime = CoarseTime[i]->at(j);
          double fineTime = FineTime[i]->at(j);
          double amp = Amplitude[i]->at(j);
          double time = coarseTime - timestamp;
          if (time > 300 && time < 500 && amp > maxAmp)
          {
            maxAmp = amp;
            index = j;
          }
        }
        if (index != -1)
        {
          RecHit rechit;
          rechit.Clear();
          rechit.setCrystalID(CrystalID);
          rechit.setFrontCenter(TVector3(Para.HitPosX(CrystalID), Para.HitPosY(CrystalID), Para.HitPosZ(CrystalID)));
          double coarseTime = CoarseTime[i]->at(index);
          double fineTime = FineTime[i]->at(index);
          double amp = Amplitude[i]->at(index);
          double time = coarseTime - timestamp - fineTime / 100;
          rechit.setEnergy(amp, Para.LGNoise(i), Para.LGPedestal(i), Para.LGMipPeak(i));
          // rechit.setEnergy(amp / correct[int(fineTime) % 1250], Para.LGNoise(i), Para.LGPedestal(i), Para.LGMipPeak(i));
          // rechit.setEnergy(amp * calibration_factor[i][int(fineTime) / 50], Para.LGNoise(i), Para.LGPedestal(i), Para.LGMipPeak(i));
          rechit.setTime(time);
          if (rechit.Energy() > 0)
          {
            HitMap[CrystalID] = rechit;
            HitID.push_back(rechit.CrystalID());
            HitTime.push_back(rechit.Time() + 817.5);
            HitEnergy.push_back(rechit.Energy());
          }
        }
      }

      Hit2Cluster m_Hit2Cluster;
      m_Hit2Cluster.Convert(HitMap, ClusterMap);
      Cluster2Shower m_Cluster2Shower;
      m_Cluster2Shower.Convert(ClusterMap, ShowerMap);

      double epsilon = 1e-10;
      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        // Calculate cluster position
        TVector3 possum;
        double etot = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          etot += it->second.Energy();
          TVector3 pos(it->second.FrontCenter());
          possum += pos * it->second.Energy();
        }
        if (etot > 0)
          possum = possum * (1 / etot);
        if (std::abs(possum.x() - (std::floor(possum.x() * 100000) / 100000)) < epsilon)
          Citer->second.setPosition(TVector3(std::floor(possum.x() * 100000) / 100000, std::floor(possum.y() * 100000) / 100000, std::floor(possum.z() * 100000) / 100000));
        else
          Citer->second.setPosition(possum);

        // Calculate cluster second moment
        double sum = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          TVector3 pos(it->second.FrontCenter());
          if (std::abs((pos - possum).Mag2()) < epsilon)
            sum += it->second.Energy() * 0;

          else
            sum += it->second.Energy() * ((pos - possum).Mag2());
        }
        if (etot > 0)
          sum /= etot;
        Citer->second.setSecondMoment(sum);
      }

      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        ClusterID.push_back(Citer->second.GetClusterID());
        ClusterSize.push_back(Citer->second.GetClusterSize());
        ClusterEnergy.push_back(Citer->second.GetClusterEnergy());
        ClusterSMoment.push_back(Citer->second.GetSecondMoment());
        ClusterPosX.push_back(Citer->second.GetPosition().x());
        ClusterPosY.push_back(Citer->second.GetPosition().y());
        ClusterPosZ.push_back(Citer->second.GetPosition().z());
      }

      for (Siter = ShowerMap.begin(); Siter != ShowerMap.end(); Siter++)
      {
        ShowerID.push_back(Siter->second.SeedID());
        ShowerE3x3.push_back(Siter->second.E3x3());
        ShowerE5x5.push_back(Siter->second.E5x5());
        ShowerEAll.push_back(Siter->second.EAll());
        ShowerPosX3x3.push_back(Siter->second.x3x3());
        ShowerPosY3x3.push_back(Siter->second.y3x3());
        ShowerPosZ3x3.push_back(Siter->second.z3x3());
        ShowerPosX5x5.push_back(Siter->second.x5x5());
        ShowerPosY5x5.push_back(Siter->second.y5x5());
        ShowerPosZ5x5.push_back(Siter->second.z5x5());
        ShowerSMoment.push_back(Siter->second.SecondMoment());
        ShowerLatMoment.push_back(Siter->second.LateralMoment());
        ShowerA20Moment.push_back(Siter->second.A20Moment());
        ShowerA42Moment.push_back(Siter->second.A42Moment());
      }
    }

    if (recmode == 2)
    {
      vector<int> SeedIDVec;
      map<int, vector<RecHit>> HitMap;
      // map<int, RecHit> HitMap;
      multimap<int, Cluster> ClusterMap;
      multimap<int, Shower> ShowerMap;

      map<int, vector<RecHit>>::iterator Hiter;
      // map<int, RecHit>::iterator Hiter;
      multimap<int, Cluster>::iterator Citer;
      multimap<int, Shower>::iterator Siter;

      SeedIDVec.clear();
      HitMap.clear();
      ClusterMap.clear();
      ShowerMap.clear();

      for (int i = 0; i < 25; i++)
      {
        int CrystalID = Hit[i]->CrystalID;
        double timestamp = TimeStamp[i];
        for (int j = 0; j < Amplitude[i]->size(); j++)
        {
          RecHit rechit;
          rechit.Clear();
          rechit.setCrystalID(CrystalID);
          rechit.setFrontCenter(TVector3(Para.HitPosX(CrystalID), Para.HitPosY(CrystalID), Para.HitPosZ(CrystalID)));
          double coarseTime = CoarseTime[i]->at(j);
          double fineTime = FineTime[i]->at(j);
          double amp = Amplitude[i]->at(j);
          double time = coarseTime - timestamp;
          rechit.setEnergy(amp, Para.HGNoise(i), Para.HGPedestal(i), Para.HGMipPeak(i));
          rechit.setTime(time);
          if (rechit.Energy() > 0)
          {
            HitMap[CrystalID].push_back(rechit);
            HitID.push_back(rechit.CrystalID());
            HitTime.push_back(rechit.Time());
            HitEnergy.push_back(rechit.Energy());
          }
        }
      }

      Hit2Cluster m_Hit2Cluster;
      m_Hit2Cluster.Convertmulti(HitMap, ClusterMap, timecut);
      Cluster2Shower m_Cluster2Shower;
      m_Cluster2Shower.Convertmulti(ClusterMap, ShowerMap);

      /*SeedFinder m_SeedFinder;
        Citer=ClusterMap.begin();
        m_SeedFinder.Seed(Citer->second,SeedIDVec);
        for(int i=0;i<SeedIDVec.size();i++) cout<<"Seed:"<<SeedIDVec.at(i)<<endl;*/

      double epsilon = 1e-10;
      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        // Calculate cluster position
        TVector3 possum;
        double etot = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          etot += it->second.Energy();
          TVector3 pos(it->second.FrontCenter());
          possum += pos * it->second.Energy();
        }
        if (etot > 0)
          possum = possum * (1 / etot);
        if (std::abs(possum.x() - (std::floor(possum.x() * 100000) / 100000)) < epsilon)
          Citer->second.setPosition(TVector3(std::floor(possum.x() * 100000) / 100000, std::floor(possum.y() * 100000) / 100000, std::floor(possum.z() * 100000) / 100000));
        else
          Citer->second.setPosition(possum);

        // Calculate cluster second moment
        double sum = 0;
        for (auto it = Citer->second.Begin(); it != Citer->second.End(); ++it)
        {
          TVector3 pos(it->second.FrontCenter());
          if (std::abs((pos - possum).Mag2()) < epsilon)
            sum += it->second.Energy() * 0;

          else
            sum += it->second.Energy() * ((pos - possum).Mag2());
        }
        if (etot > 0)
          sum /= etot;
        Citer->second.setSecondMoment(sum);
      }

      for (Citer = ClusterMap.begin(); Citer != ClusterMap.end(); Citer++)
      {
        ClusterID.push_back(Citer->second.GetClusterID());
        ClusterSize.push_back(Citer->second.GetClusterSize());
        ClusterEnergy.push_back(Citer->second.GetClusterEnergy());
        ClusterSMoment.push_back(Citer->second.GetSecondMoment());
        ClusterPosX.push_back(Citer->second.GetPosition().x());
        ClusterPosY.push_back(Citer->second.GetPosition().y());
        ClusterPosZ.push_back(Citer->second.GetPosition().z());
      }

      for (Siter = ShowerMap.begin(); Siter != ShowerMap.end(); Siter++)
      {
        ShowerID.push_back(Siter->second.SeedID());
        ShowerE3x3.push_back(Siter->second.E3x3());
        ShowerE5x5.push_back(Siter->second.E5x5());
        ShowerEAll.push_back(Siter->second.EAll());
        ShowerPosX3x3.push_back(Siter->second.x3x3());
        ShowerPosY3x3.push_back(Siter->second.y3x3());
        ShowerPosZ3x3.push_back(Siter->second.z3x3());
        ShowerPosX5x5.push_back(Siter->second.x5x5());
        ShowerPosY5x5.push_back(Siter->second.y5x5());
        ShowerPosZ5x5.push_back(Siter->second.z5x5());
        ShowerSMoment.push_back(Siter->second.SecondMoment());
        ShowerLatMoment.push_back(Siter->second.LateralMoment());
        ShowerA20Moment.push_back(Siter->second.A20Moment());
        ShowerA42Moment.push_back(Siter->second.A42Moment());
      }
    }

    rec_data->Fill();
  }

  f->cd();
  rec_data->Write();
  f->Close();

  for (int i = 0; i < 25; i++)
  {
    delete Hit[i];
    delete CoarseTime[i];
    delete FineTime[i];
    delete Amplitude[i];
  }

  return 0;
}
