/*
 * File: cunhj.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cunhj.h"
#include "sqrt1.h"
#include "cmlri.h"
#include "log1.h"
#include "Calc_pressao_estacionario_rtwutil.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int ipmtr
 *                double tol
 *                creal_T *phi
 *                creal_T *arg
 *                creal_T *zeta1
 *                creal_T *zeta2
 *                creal_T *asum
 *                creal_T *bsum
 * Return Type  : void
 */
void b_cunhj(const creal_T z, double fnu, int ipmtr, double tol, creal_T *phi,
             creal_T *arg, creal_T *zeta1, creal_T *zeta2, creal_T *asum,
             creal_T *bsum)
{
  double rfnu;
  creal_T up[14];
  double ac;
  creal_T zb;
  double rfnu2;
  double fn23;
  double rfn13_re;
  creal_T w2;
  int k;
  creal_T w;
  creal_T p[30];
  double ap[30];
  double pp;
  creal_T suma;
  int l1;
  boolean_T exitg1;
  double brm;
  creal_T zc;
  double bim;
  double d;
  static const double GAMA[30] = { 0.6299605249474366, 0.25198420997897464,
    0.15479030041565583, 0.11071306241615901, 0.085730939552739485,
    0.069716131695868433, 0.058608567189371359, 0.050469887353631067,
    0.044260058068915482, 0.039372066154350994, 0.035428319592445537,
    0.032181885750209825, 0.029464624079115768, 0.027158167711293448,
    0.025176827297386177, 0.023457075530607888, 0.021950839013490719,
    0.020621082823564625, 0.019438824089788084, 0.018381063380068317,
    0.017429321323196318, 0.016568583778661234, 0.015786528598791844,
    0.015072950149409559, 0.014419325083995464, 0.013818480573534178,
    0.013264337899427657, 0.012751712197049864, 0.012276154531876277,
    0.011833826239848241 };

  double zth_re;
  double zth_im;
  static const double BETA[210] = { 0.017998872141355329, 0.0055996491106438812,
    0.0028850140223113277, 0.0018009660676105393, 0.0012475311058919921,
    0.00092287887657293828, 0.00071443042172728737, 0.00057178728178970488,
    0.00046943100760648155, 0.00039323283546291665, 0.00033481888931829768,
    0.00028895214849575154, 0.0002522116155495733, 0.00022228058079888332,
    0.00019754183803306251, 0.00017683685501971802, 0.00015931689966182109,
    0.00014434793019733397, 0.00013144806811996539, 0.00012024544494930288,
    0.00011044914450459939, 0.00010182877074056726, 9.4199822420423752E-5,
    8.7413054575383449E-5, 8.1346626216280142E-5, 7.590022696462193E-5,
    7.0990630063415351E-5, 6.6548287484246817E-5, 6.25146958969275E-5,
    5.8840339442625178E-5, -0.0014928295321342917, -0.00087820470954638936,
    -0.00050291654957203464, -0.000294822138512746, -0.00017546399697078284,
    -0.00010400855046081644, -5.961419530464579E-5, -3.1203892907609836E-5,
    -1.2608973598023005E-5, -2.4289260857573037E-7, 8.059961654142736E-6,
    1.3650700926214739E-5, 1.7396412547292627E-5, 1.9867297884213378E-5,
    2.1446326379082263E-5, 2.2395465923245652E-5, 2.2896778381471263E-5,
    2.3078538981117782E-5, 2.3032197608090914E-5, 2.2823607372034874E-5,
    2.2500588110529241E-5, 2.2098101536199144E-5, 2.164184274481039E-5,
    2.1150764925622083E-5, 2.0638874978217072E-5, 2.0116524199708165E-5,
    1.9591345014117925E-5, 1.9068936791043675E-5, 1.8553371964163667E-5,
    1.8047572225967421E-5, 0.0005522130767212928, 0.00044793258155238465,
    0.00027952065399202059, 0.0001524681561984466, 6.932711056570436E-5,
    1.7625868306999139E-5, -1.3574499634326914E-5, -3.1797241335042717E-5,
    -4.188618616966934E-5, -4.6900488937914104E-5, -4.8766544741378735E-5,
    -4.8701003118673505E-5, -4.7475562089008661E-5, -4.5581305813862843E-5,
    -4.33309644511266E-5, -4.0923019315775034E-5, -3.848226386032213E-5,
    -3.6085716753541052E-5, -3.3779330612336739E-5, -3.1588856077210961E-5,
    -2.9526956175080731E-5, -2.7597891482833575E-5, -2.5800617466688372E-5,
    -2.4130835676128019E-5, -2.2582350951834605E-5, -2.1147965676891298E-5,
    -1.9820063888529493E-5, -1.8590987080106508E-5, -1.7453269984421023E-5,
    -1.63997823854498E-5, -0.0004746177965599598, -0.0004778645671473215,
    -0.00032039022806703763, -0.00016110501611996228, -4.2577810128543523E-5,
    3.4457129429496748E-5, 7.97092684075675E-5, 0.00010313823670827221,
    0.00011246677526220416, 0.00011310364210848139, 0.00010865163484877427,
    0.00010143795159766197, 9.29298396593364E-5, 8.4029313301609E-5,
    7.52727991349134E-5, 6.696325219757309E-5, 5.925645473231947E-5,
    5.2216930882697554E-5, 4.5853948516536063E-5, 4.0144551389148682E-5,
    3.5048173003132809E-5, 3.0515799503434667E-5, 2.6495611995051603E-5,
    2.2936363369099816E-5, 1.9789305666402162E-5, 1.7009198463641262E-5,
    1.45547428261524E-5, 1.2388664099587841E-5, 1.0477587607658323E-5,
    8.7917995497847932E-6, 0.00073646581057257841, 0.000872790805146194,
    0.00062261486257313506, 0.00028599815419430417, 3.8473767287936606E-6,
    -0.00018790600363697156, -0.00029760364659455455, -0.00034599812683265633,
    -0.00035338247091603773, -0.00033571563577504876, -0.00030432112478903981,
    -0.00026672272304761283, -0.00022765421412281953, -0.00018992261185456235,
    -0.00015505891859909386, -0.00012377824076187363, -9.6292614771764412E-5,
    -7.2517832771442527E-5, -5.2207002889563382E-5, -3.5034775051190054E-5,
    -2.0648976103555174E-5, -8.7010609684976711E-6, 1.136986866751003E-6,
    9.1642647412277879E-6, 1.5647778542887261E-5, 2.0822362948246685E-5,
    2.4892338100459516E-5, 2.8034050957414632E-5, 3.0398777462986191E-5,
    3.2115673140670063E-5, -0.0018018219196388571, -0.0024340296293804253,
    -0.001834226635498568, -0.00076220459635400974, 0.00023907947525692722,
    0.00094926611717688109, 0.0013446744970154036, 0.0014845749525944918,
    0.0014473233983061759, 0.0013026826128565718, 0.0011035159737564268,
    0.00088604744041979172, 0.00067307320816566542, 0.00047760387285658237,
    0.00030599192635878935, 0.00016031569459472162, 4.0074955527061327E-5,
    -5.6660746163525162E-5, -0.00013250618677298264, -0.00019029618798961406,
    -0.00023281145037693741, -0.00026262881146466884, -0.00028205046986759866,
    -0.00029308156319286116, -0.00029743596217631662, -0.00029655733423934809,
    -0.00029164736331209088, -0.00028369620383773418, -0.00027351231709567335,
    -0.00026175015580676858, 0.0063858589121205088, 0.00962374215806378,
    0.0076187806120700105, 0.0028321905554562804, -0.0020984135201272008,
    -0.0057382676421662646, -0.0077080424449541465, -0.0082101169226484437,
    -0.0076582452034690543, -0.006472097293910452, -0.0049913241200496648,
    -0.0034561228971313326, -0.0020178558001417079, -0.00075943068678196145,
    0.00028417363152385912, 0.0011089166758633741, 0.0017290149387272878,
    0.0021681259080268472, 0.0024535771049453972, 0.0026128182105833488,
    0.0026714103965627691, 0.0026520307339598045, 0.0025741165287728731,
    0.0024538912623609443, 0.002304600580717955, 0.0021368483768671267,
    0.0019589652847887091, 0.0017773700867945441, 0.0015969028076583906,
    0.0014211197566443854 };

  int l2;
  int ias;
  int ibs;
  int is;
  double rtzta_re;
  double rtzta_im;
  int ks;
  boolean_T exitg2;
  static const double ALFA[180] = { -0.0044444444444444444,
    -0.000922077922077922, -8.8489288489288488E-5, 0.00016592768783244973,
    0.00024669137274179289, 0.00026599558934625478, 0.00026182429706150096,
    0.00024873043734465562, 0.00023272104008323209, 0.00021636248571236508,
    0.00020073885876275234, 0.00018626763663754517, 0.0001730607759178765,
    0.00016109170592901574, 0.00015027477416090814, 0.00014050349739126979,
    0.0001316688165459228, 0.00012366744559825325, 0.00011640527147473791,
    0.00010979829837271337, 0.00010377241042299283, 9.8262607836936344E-5,
    9.321205172495032E-5, 8.857108524787117E-5, 8.4296310571570029E-5,
    8.0349754840779115E-5, 7.6698134535920737E-5, 7.3312215748177779E-5,
    7.0166262516314139E-5, 6.7237563379016026E-5, 0.000693735541354589,
    0.00023224174518292166, -1.419862735566912E-5, -0.00011644493167204864,
    -0.00015080355805304876, -0.00015512192491809622, -0.00014680975664646556,
    -0.00013381550386749137, -0.00011974497568425405, -0.00010618431920797402,
    -9.3769954989119444E-5, -8.2692304558819327E-5, -7.2937434815522126E-5,
    -6.4404235772101633E-5, -5.69611566009369E-5, -5.0473104430356164E-5,
    -4.4813486800888282E-5, -3.9868872771759884E-5, -3.5540053297204251E-5,
    -3.174142566090225E-5, -2.839967939041748E-5, -2.5452272063487058E-5,
    -2.2845929716472455E-5, -2.0535275310648061E-5, -1.848162176276661E-5,
    -1.6651933002139381E-5, -1.5017941298011949E-5, -1.3555403137904052E-5,
    -1.2243474647385812E-5, -1.1064188481130817E-5, -0.00035421197145774384,
    -0.00015616126394515941, 3.0446550359493642E-5, 0.00013019865577324269,
    0.00016747110669971228, 0.00017022258768359256, 0.00015650142760859472,
    0.00013633917097744512, 0.00011488669202982512, 9.4586909303468817E-5,
    7.6449841925089825E-5, 6.0757033496519734E-5, 4.7439429929050881E-5,
    3.6275751200534429E-5, 2.6993971497922491E-5, 1.9321093824793926E-5,
    1.3005667479396321E-5, 7.8262086674449658E-6, 3.5925748581935159E-6,
    1.4404004981425182E-7, -2.6539676969793912E-6, -4.9134686709848593E-6,
    -6.7273929609124832E-6, -8.17269379678658E-6, -9.313047150935612E-6,
    -1.0201141879801643E-5, -1.0880596251059288E-5, -1.1387548150960355E-5,
    -1.1751967567455642E-5, -1.1998736487094414E-5, 0.00037819419920177291,
    0.00020247195276181616, -6.3793850631886236E-5, -0.0002385982306030059,
    -0.00031091625602736159, -0.00031368011524757634, -0.00027895027379132341,
    -0.00022856408261914138, -0.00017524528034084676, -0.00012554406306069035,
    -8.2298287282020835E-5, -4.6286073058811649E-5, -1.7233430236696227E-5,
    5.6069048230460226E-6, 2.313954431482868E-5, 3.6264274585679393E-5,
    4.5800612449018877E-5, 5.2459529495911405E-5, 5.6839620854581527E-5,
    5.9434982039310406E-5, 6.0647852757842175E-5, 6.0802390778843649E-5,
    6.0157789453946036E-5, 5.8919965734469847E-5, 5.72515823777593E-5,
    5.5280437558585257E-5, 5.3106377380288019E-5, 5.080693020123257E-5,
    4.8441864762009484E-5, 4.6056858160747536E-5, -0.00069114139728829421,
    -0.00042997663305887192, 0.000183067735980039, 0.00066008814754201417,
    0.0008759649699511859, 0.00087733523595823551, 0.00074936958537899067,
    0.000563832329756981, 0.00036805931997144317, 0.00018846453551445559,
    3.7066305766490415E-5, -8.28520220232137E-5, -0.000172751952869173,
    -0.00023631487360587297, -0.00027796615069490668, -0.00030207951415545694,
    -0.00031259471264382012, -0.00031287255875806717, -0.0003056780384663244,
    -0.00029322647061455731, -0.0002772556555829348, -0.00025910392846703172,
    -0.00023978401439648034, -0.00022004826004542284, -0.00020044391109497149,
    -0.00018135869221097068, -0.00016305767447865748, -0.00014571267217520584,
    -0.0001294254219839246, -0.00011424569194244596, 0.0019282196424877589,
    0.0013559257630202223, -0.000717858090421303, -0.0025808480257527035,
    -0.0034927113082616847, -0.0034698629934096061, -0.0028228523335131019,
    -0.0018810307640489134, -0.00088953171838394764, 3.8791210263103525E-6,
    0.00072868854011969139, 0.0012656637305345775, 0.0016251815837267443,
    0.0018320315321637317, 0.0019158838899052792, 0.0019058884675554615,
    0.0018279898242182574, 0.0017038950642112153, 0.0015509712717109768,
    0.0013826142185227616, 0.0012088142423006478, 0.0010367653263834496,
    0.00087143791806861914, 0.000716080155297701, 0.00057263700255812935,
    0.00044208981946580229, 0.00032472494850309055, 0.00022034204273024659,
    0.00012841289840135388, 4.8200592455209545E-5 };

  double tfn_re;
  double tfn_im;
  double rzth_re;
  double rzth_im;
  int kp1;
  int l;
  creal_T cr[14];
  creal_T dr[14];
  static const double C[105] = { 1.0, -0.20833333333333334, 0.125,
    0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
    1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
    -11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
    -28.212072558200244, 84.636217674600729, -91.818241543240021,
    42.534998745388457, -7.3687943594796321, 0.22710800170898438,
    212.57013003921713, -765.25246814118168, 1059.9904525279999,
    -699.57962737613252, 218.19051174421159, -26.491430486951554,
    0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
    -13586.550006434138, 11655.393336864534, -5305.646978613403,
    1200.9029132163525, -108.09091978839466, 1.7277275025844574,
    20204.291330966149, -96980.598388637518, 192547.00123253153,
    -203400.17728041555, 122200.46498301746, -41192.65496889755,
    7109.5143024893641, -493.915304773088, 6.074042001273483,
    -242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
    3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
    -331645.17248456361, 45218.768981362729, -2499.8304818112097,
    24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
    5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
    -3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
    308186.40461266239, -13886.08975371704, 110.01714026924674,
    -4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
    1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
    -4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
    2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
    8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
    -3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
    1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
    -2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
    3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
    -3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
    1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
    -8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
    1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175 };

  static const double BR[14] = { 1.0, -0.14583333333333334,
    -0.098741319444444448, -0.14331205391589505, -0.31722720267841353,
    -0.94242914795712029, -3.5112030408263544, -15.727263620368046,
    -82.281439097185938, -492.3553705236705, -3316.2185685479726,
    -24827.674245208589, -204526.58731512979, -1.83844491706821E+6 };

  static const double AR[14] = { 1.0, 0.10416666666666667, 0.083550347222222224,
    0.12822657455632716, 0.29184902646414046, 0.88162726744375763,
    3.3214082818627677, 14.995762986862555, 78.923013011586519,
    474.45153886826432, 3207.490090890662, 24086.549640874004, 198923.1191695098,
    1.7919020077753437E+6 };

  rfnu = 1.0 / fnu;
  memset(&up[0], 0, 14U * sizeof(creal_T));
  ac = fnu * 2.2250738585072014E-305;
  asum->re = 0.0;
  asum->im = 0.0;
  bsum->re = 0.0;
  bsum->im = 0.0;
  if ((fabs(z.re) <= ac) && (fabs(z.im) <= ac)) {
    zeta1->re = 1402.9773265065639 + fnu;
    zeta1->im = 0.0;
    zeta2->re = fnu;
    zeta2->im = 0.0;
    phi->re = 1.0;
    phi->im = 0.0;
    arg->re = 1.0;
    arg->im = 0.0;
  } else {
    zb.re = rfnu * z.re;
    zb.im = rfnu * z.im;
    rfnu2 = rfnu * rfnu;
    ac = rt_powd_snf(fnu, 0.33333333333333331);
    fn23 = ac * ac;
    rfn13_re = 1.0 / ac;
    w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
    w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
    ac = rt_hypotd_snf(w2.re, w2.im);
    if (ac > 0.25) {
      w = w2;
      b_sqrt(&w);
      ac = w.re;
      pp = w.im;
      if (w.re < 0.0) {
        ac = 0.0;
      }

      if (w.im < 0.0) {
        pp = 0.0;
      }

      w.re = ac;
      w.im = pp;
      if (zb.im == 0.0) {
        if (pp == 0.0) {
          zc.re = (1.0 + ac) / zb.re;
          zc.im = 0.0;
        } else if (1.0 + ac == 0.0) {
          zc.re = 0.0;
          zc.im = pp / zb.re;
        } else {
          zc.re = (1.0 + ac) / zb.re;
          zc.im = pp / zb.re;
        }
      } else if (zb.re == 0.0) {
        if (1.0 + ac == 0.0) {
          zc.re = pp / zb.im;
          zc.im = 0.0;
        } else if (pp == 0.0) {
          zc.re = 0.0;
          zc.im = -((1.0 + ac) / zb.im);
        } else {
          zc.re = pp / zb.im;
          zc.im = -((1.0 + ac) / zb.im);
        }
      } else {
        brm = fabs(zb.re);
        bim = fabs(zb.im);
        if (brm > bim) {
          bim = zb.im / zb.re;
          d = zb.re + bim * zb.im;
          zc.re = ((1.0 + ac) + bim * pp) / d;
          zc.im = (pp - bim * (1.0 + ac)) / d;
        } else if (bim == brm) {
          if (zb.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (zb.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          zc.re = ((1.0 + ac) * bim + pp * d) / brm;
          zc.im = (pp * bim - (1.0 + ac) * d) / brm;
        } else {
          bim = zb.re / zb.im;
          d = zb.im + bim * zb.re;
          zc.re = (bim * (1.0 + ac) + pp) / d;
          zc.im = (bim * pp - (1.0 + ac)) / d;
        }
      }

      b_log(&zc);
      ac = zc.re;
      pp = zc.im;
      if (zc.re < 0.0) {
        ac = 0.0;
      }

      if (zc.im < 0.0) {
        pp = 0.0;
      } else {
        if (zc.im > 1.5707963267948966) {
          pp = 1.5707963267948966;
        }
      }

      zth_re = 1.5 * (ac - w.re);
      zth_im = 1.5 * (pp - w.im);
      zeta1->re = fnu * ac;
      zeta1->im = fnu * pp;
      zeta2->re = fnu * w.re;
      zeta2->im = fnu * w.im;
      if ((zth_re >= 0.0) && (zth_im < 0.0)) {
        ac = 4.71238898038469;
      } else if (zth_re != 0.0) {
        ac = atan(zth_im / zth_re);
        if (zth_re < 0.0) {
          ac += 3.1415926535897931;
        }
      } else {
        ac = 1.5707963267948966;
      }

      pp = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
      ac *= 0.66666666666666663;
      zb.re = pp * cos(ac);
      zb.im = pp * sin(ac);
      if (zb.im < 0.0) {
        zb.im = 0.0;
      }

      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      if (zb.im == 0.0) {
        if (zth_im == 0.0) {
          rtzta_re = zth_re / zb.re;
          rtzta_im = 0.0;
        } else if (zth_re == 0.0) {
          rtzta_re = 0.0;
          rtzta_im = zth_im / zb.re;
        } else {
          rtzta_re = zth_re / zb.re;
          rtzta_im = zth_im / zb.re;
        }
      } else if (zb.re == 0.0) {
        if (zth_re == 0.0) {
          rtzta_re = zth_im / zb.im;
          rtzta_im = 0.0;
        } else if (zth_im == 0.0) {
          rtzta_re = 0.0;
          rtzta_im = -(zth_re / zb.im);
        } else {
          rtzta_re = zth_im / zb.im;
          rtzta_im = -(zth_re / zb.im);
        }
      } else {
        brm = fabs(zb.re);
        bim = zb.im;
        if (brm > bim) {
          bim = zb.im / zb.re;
          d = zb.re + bim * zb.im;
          rtzta_re = (zth_re + bim * zth_im) / d;
          rtzta_im = (zth_im - bim * zth_re) / d;
        } else if (bim == brm) {
          if (zb.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (zb.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          rtzta_re = (zth_re * bim + zth_im * d) / brm;
          rtzta_im = (zth_im * bim - zth_re * d) / brm;
        } else {
          bim = zb.re / zb.im;
          d = zb.im + bim * zb.re;
          rtzta_re = (bim * zth_re + zth_im) / d;
          rtzta_im = (bim * zth_im - zth_re) / d;
        }
      }

      if (w.im == 0.0) {
        if (rtzta_im == 0.0) {
          suma.re = rtzta_re / w.re;
          suma.im = 0.0;
        } else if (rtzta_re == 0.0) {
          suma.re = 0.0;
          suma.im = rtzta_im / w.re;
        } else {
          suma.re = rtzta_re / w.re;
          suma.im = rtzta_im / w.re;
        }
      } else if (w.re == 0.0) {
        if (rtzta_re == 0.0) {
          suma.re = rtzta_im / w.im;
          suma.im = 0.0;
        } else if (rtzta_im == 0.0) {
          suma.re = 0.0;
          suma.im = -(rtzta_re / w.im);
        } else {
          suma.re = rtzta_im / w.im;
          suma.im = -(rtzta_re / w.im);
        }
      } else {
        brm = fabs(w.re);
        bim = fabs(w.im);
        if (brm > bim) {
          bim = w.im / w.re;
          d = w.re + bim * w.im;
          suma.re = (rtzta_re + bim * rtzta_im) / d;
          suma.im = (rtzta_im - bim * rtzta_re) / d;
        } else if (bim == brm) {
          if (w.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (w.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          suma.re = (rtzta_re * bim + rtzta_im * d) / brm;
          suma.im = (rtzta_im * bim - rtzta_re * d) / brm;
        } else {
          bim = w.re / w.im;
          d = w.im + bim * w.re;
          suma.re = (bim * rtzta_re + rtzta_im) / d;
          suma.im = (bim * rtzta_im - rtzta_re) / d;
        }
      }

      zb.re = suma.re + suma.re;
      zb.im = suma.im + suma.im;
      b_sqrt(&zb);
      phi->re = zb.re * rfn13_re - zb.im * 0.0;
      phi->im = zb.re * 0.0 + zb.im * rfn13_re;
      if (ipmtr == 1) {
      } else {
        if (w.im == 0.0) {
          tfn_re = rfnu / w.re;
          tfn_im = 0.0;
        } else if (w.re == 0.0) {
          if (rfnu == 0.0) {
            tfn_re = 0.0 / w.im;
            tfn_im = 0.0;
          } else {
            tfn_re = 0.0;
            tfn_im = -(rfnu / w.im);
          }
        } else {
          brm = fabs(w.re);
          bim = fabs(w.im);
          if (brm > bim) {
            bim = w.im / w.re;
            d = w.re + bim * w.im;
            tfn_re = (rfnu + bim * 0.0) / d;
            tfn_im = (0.0 - bim * rfnu) / d;
          } else if (bim == brm) {
            if (w.re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (w.im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            tfn_re = (rfnu * bim + 0.0 * d) / brm;
            tfn_im = (0.0 * bim - rfnu * d) / brm;
          } else {
            bim = w.re / w.im;
            d = w.im + bim * w.re;
            tfn_re = bim * rfnu / d;
            tfn_im = (bim * 0.0 - rfnu) / d;
          }
        }

        if (zth_im == 0.0) {
          rzth_re = rfnu / zth_re;
          rzth_im = 0.0;
        } else if (zth_re == 0.0) {
          if (rfnu == 0.0) {
            rzth_re = 0.0 / zth_im;
            rzth_im = 0.0;
          } else {
            rzth_re = 0.0;
            rzth_im = -(rfnu / zth_im);
          }
        } else {
          brm = fabs(zth_re);
          bim = fabs(zth_im);
          if (brm > bim) {
            bim = zth_im / zth_re;
            d = zth_re + bim * zth_im;
            rzth_re = (rfnu + bim * 0.0) / d;
            rzth_im = (0.0 - bim * rfnu) / d;
          } else if (bim == brm) {
            if (zth_re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (zth_im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            rzth_re = (rfnu * bim + 0.0 * d) / brm;
            rzth_im = (0.0 * bim - rfnu * d) / brm;
          } else {
            bim = zth_re / zth_im;
            d = zth_im + bim * zth_re;
            rzth_re = bim * rfnu / d;
            rzth_im = (bim * 0.0 - rfnu) / d;
          }
        }

        zc.re = 0.10416666666666667 * rzth_re;
        zc.im = 0.10416666666666667 * rzth_im;
        if (w2.im == 0.0) {
          w.re = 1.0 / w2.re;
          w.im = 0.0;
        } else if (w2.re == 0.0) {
          w.re = 0.0;
          w.im = -(1.0 / w2.im);
        } else {
          brm = fabs(w2.re);
          bim = fabs(w2.im);
          if (brm > bim) {
            bim = w2.im / w2.re;
            d = w2.re + bim * w2.im;
            w.re = (1.0 + bim * 0.0) / d;
            w.im = (0.0 - bim) / d;
          } else if (bim == brm) {
            if (w2.re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (w2.im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            w.re = (bim + 0.0 * d) / brm;
            w.im = (0.0 * bim - d) / brm;
          } else {
            bim = w2.re / w2.im;
            d = w2.im + bim * w2.re;
            w.re = bim / d;
            w.im = (bim * 0.0 - 1.0) / d;
          }
        }

        ac = w.re * -0.20833333333333334 + 0.125;
        pp = w.im * -0.20833333333333334;
        up[1].re = ac * tfn_re - pp * tfn_im;
        up[1].im = ac * tfn_im + pp * tfn_re;
        bsum->re = up[1].re + zc.re;
        bsum->im = up[1].im + zc.im;
        if (rfnu >= tol) {
          zth_re = rzth_re;
          zth_im = rzth_im;
          w2.re = tfn_re;
          w2.im = tfn_im;
          up[0].re = 1.0;
          up[0].im = 0.0;
          pp = 1.0;
          fn23 = tol * (fabs(bsum->re) + fabs(bsum->im));
          ks = -1;
          kp1 = 2;
          l = 2;
          ias = 0;
          ibs = 0;
          for (l1 = 0; l1 < 14; l1++) {
            cr[l1].re = 0.0;
            cr[l1].im = 0.0;
            dr[l1].re = 0.0;
            dr[l1].im = 0.0;
          }

          is = 2;
          exitg1 = false;
          while ((!exitg1) && (is < 13)) {
            for (k = is; k <= is + 1; k++) {
              ks++;
              kp1++;
              l++;
              suma.re = C[l];
              suma.im = 0.0;
              for (l1 = 2; l1 <= kp1; l1++) {
                l++;
                ac = suma.re;
                suma.re = (suma.re * w.re - suma.im * w.im) + C[l];
                suma.im = ac * w.im + suma.im * w.re;
              }

              ac = w2.re;
              w2.re = w2.re * tfn_re - w2.im * tfn_im;
              w2.im = ac * tfn_im + w2.im * tfn_re;
              up[kp1 - 1].re = w2.re * suma.re - w2.im * suma.im;
              up[kp1 - 1].im = w2.re * suma.im + w2.im * suma.re;
              cr[ks].re = BR[ks + 1] * zth_re;
              cr[ks].im = BR[ks + 1] * zth_im;
              ac = zth_re;
              zth_re = zth_re * rzth_re - zth_im * rzth_im;
              zth_im = ac * rzth_im + zth_im * rzth_re;
              dr[ks].re = AR[ks + 2] * zth_re;
              dr[ks].im = AR[ks + 2] * zth_im;
            }

            pp *= rfnu2;
            if (ias != 1) {
              suma = up[is];
              l1 = is;
              for (l2 = 1; l2 <= is; l2++) {
                l1--;
                suma.re += cr[l2 - 1].re * up[l1].re - cr[l2 - 1].im * up[l1].im;
                suma.im += cr[l2 - 1].re * up[l1].im + cr[l2 - 1].im * up[l1].re;
              }

              asum->re += suma.re;
              asum->im += suma.im;
              if ((pp < tol) && (fabs(asum->re) + fabs(asum->im) < tol)) {
                ias = 1;
              }
            }

            if (ibs != 1) {
              zb.re = up[is + 1].re + (up[is].re * zc.re - up[is].im * zc.im);
              zb.im = up[is + 1].im + (up[is].re * zc.im + up[is].im * zc.re);
              l1 = is;
              for (l2 = 1; l2 <= is; l2++) {
                l1--;
                zb.re += dr[l2 - 1].re * up[l1].re - dr[l2 - 1].im * up[l1].im;
                zb.im += dr[l2 - 1].re * up[l1].im + dr[l2 - 1].im * up[l1].re;
              }

              bsum->re += zb.re;
              bsum->im += zb.im;
              if ((pp < fn23) && (fabs(bsum->re) + fabs(bsum->im) < tol)) {
                ibs = 1;
              }
            }

            if ((ias == 1) && (ibs == 1)) {
              exitg1 = true;
            } else {
              is += 2;
            }
          }
        }

        asum->re++;
        bsum->re = -bsum->re;
        bsum->im = -bsum->im;
        ac = bsum->re;
        pp = bsum->im;
        bsum->re = bsum->re * rfn13_re - bsum->im * 0.0;
        bsum->im = ac * 0.0 + bsum->im * rfn13_re;
        fn23 = ac * rfn13_re - pp * 0.0;
        pp = ac * 0.0 + pp * rfn13_re;
        if (rtzta_im == 0.0) {
          if (bsum->im == 0.0) {
            bsum->re /= rtzta_re;
            bsum->im = 0.0;
          } else if (bsum->re == 0.0) {
            bsum->re = 0.0;
            bsum->im /= rtzta_re;
          } else {
            bsum->re /= rtzta_re;
            bsum->im /= rtzta_re;
          }
        } else if (rtzta_re == 0.0) {
          if (bsum->re == 0.0) {
            bsum->re = bsum->im / rtzta_im;
            bsum->im = 0.0;
          } else if (bsum->im == 0.0) {
            bsum->re = 0.0;
            bsum->im = -(fn23 / rtzta_im);
          } else {
            bsum->re = bsum->im / rtzta_im;
            bsum->im = -(fn23 / rtzta_im);
          }
        } else {
          brm = fabs(rtzta_re);
          bim = fabs(rtzta_im);
          if (brm > bim) {
            bim = rtzta_im / rtzta_re;
            d = rtzta_re + bim * rtzta_im;
            bsum->re = (fn23 + bim * pp) / d;
            bsum->im = (pp - bim * fn23) / d;
          } else if (bim == brm) {
            if (rtzta_re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (rtzta_im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            bsum->re = (fn23 * bim + pp * d) / brm;
            bsum->im = (pp * bim - fn23 * d) / brm;
          } else {
            bim = rtzta_re / rtzta_im;
            d = rtzta_im + bim * rtzta_re;
            bsum->re = (bim * fn23 + pp) / d;
            bsum->im = (bim * pp - fn23) / d;
          }
        }
      }
    } else {
      k = 1;
      memset(&p[0], 0, 30U * sizeof(creal_T));
      memset(&ap[0], 0, 30U * sizeof(double));
      p[0].re = 1.0;
      p[0].im = 0.0;
      suma.re = 0.6299605249474366;
      suma.im = 0.0;
      ap[0] = 1.0;
      if (ac >= tol) {
        k = 30;
        l1 = 1;
        exitg1 = false;
        while ((!exitg1) && (l1 + 1 < 31)) {
          pp = p[l1 - 1].re;
          p[l1].re = p[l1 - 1].re * w2.re - p[l1 - 1].im * w2.im;
          p[l1].im = pp * w2.im + p[l1 - 1].im * w2.re;
          suma.re += p[l1].re * GAMA[l1];
          suma.im += p[l1].im * GAMA[l1];
          ap[l1] = ap[l1 - 1] * ac;
          if (ap[l1] < tol) {
            k = l1 + 1;
            exitg1 = true;
          } else {
            l1++;
          }
        }
      }

      zb.re = w2.re * suma.re - w2.im * suma.im;
      zb.im = w2.re * suma.im + w2.im * suma.re;
      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      b_sqrt(&suma);
      b_sqrt(&w2);
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      ac = (zb.re * suma.re - zb.im * suma.im) * 0.66666666666666663 + 1.0;
      pp = (zb.re * suma.im + zb.im * suma.re) * 0.66666666666666663;
      zeta1->re = zeta2->re * ac - zeta2->im * pp;
      zeta1->im = zeta2->re * pp + zeta2->im * ac;
      suma.re += suma.re;
      suma.im += suma.im;
      b_sqrt(&suma);
      phi->re = suma.re * rfn13_re - suma.im * 0.0;
      phi->im = suma.re * 0.0 + suma.im * rfn13_re;
      if (ipmtr == 1) {
      } else {
        zb.re = 0.0;
        zb.im = 0.0;
        for (l1 = 0; l1 + 1 <= k; l1++) {
          zb.re += p[l1].re * BETA[l1];
          zb.im += p[l1].im * BETA[l1];
        }

        *bsum = zb;
        l1 = 0;
        l2 = 30;
        fn23 = tol * rt_hypotd_snf(zb.re, zb.im);
        ac = tol;
        pp = 1.0;
        ias = 0;
        ibs = 0;
        if (rfnu2 >= tol) {
          is = 2;
          exitg1 = false;
          while ((!exitg1) && (is < 8)) {
            ac /= rfnu2;
            pp *= rfnu2;
            if (ias != 1) {
              suma.re = 0.0;
              suma.im = 0.0;
              ks = 0;
              exitg2 = false;
              while ((!exitg2) && (ks + 1 <= k)) {
                suma.re += p[ks].re * ALFA[l1 + ks];
                suma.im += p[ks].im * ALFA[l1 + ks];
                if (ap[ks] < ac) {
                  exitg2 = true;
                } else {
                  ks++;
                }
              }

              asum->re += suma.re * pp;
              asum->im += suma.im * pp;
              if (pp < tol) {
                ias = 1;
              }
            }

            if (ibs != 1) {
              zb.re = 0.0;
              zb.im = 0.0;
              ks = 0;
              exitg2 = false;
              while ((!exitg2) && (ks + 1 <= k)) {
                zb.re += p[ks].re * BETA[l2 + ks];
                zb.im += p[ks].im * BETA[l2 + ks];
                if (ap[ks] < ac) {
                  exitg2 = true;
                } else {
                  ks++;
                }
              }

              bsum->re += zb.re * pp;
              bsum->im += zb.im * pp;
              if (pp < fn23) {
                ibs = 1;
              }
            }

            if ((ias == 1) && (ibs == 1)) {
              exitg1 = true;
            } else {
              l1 += 30;
              l2 += 30;
              is++;
            }
          }
        }

        asum->re++;
        bsum->re *= rfnu * rfn13_re;
        bsum->im *= rfnu * rfn13_re;
      }
    }
  }
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                double tol
 *                creal_T *phi
 *                creal_T *arg
 *                creal_T *zeta1
 *                creal_T *zeta2
 * Return Type  : void
 */
void cunhj(const creal_T z, double fnu, double tol, creal_T *phi, creal_T *arg,
           creal_T *zeta1, creal_T *zeta2)
{
  double ac;
  creal_T zb;
  double fn23;
  double rfn13_re;
  creal_T w2;
  creal_T p[30];
  double ap[30];
  double tmpr;
  double tmpi;
  creal_T suma;
  int i;
  boolean_T exitg1;
  double brm;
  creal_T zc;
  double pp;
  static const double GAMA[30] = { 0.6299605249474366, 0.25198420997897464,
    0.15479030041565583, 0.11071306241615901, 0.085730939552739485,
    0.069716131695868433, 0.058608567189371359, 0.050469887353631067,
    0.044260058068915482, 0.039372066154350994, 0.035428319592445537,
    0.032181885750209825, 0.029464624079115768, 0.027158167711293448,
    0.025176827297386177, 0.023457075530607888, 0.021950839013490719,
    0.020621082823564625, 0.019438824089788084, 0.018381063380068317,
    0.017429321323196318, 0.016568583778661234, 0.015786528598791844,
    0.015072950149409559, 0.014419325083995464, 0.013818480573534178,
    0.013264337899427657, 0.012751712197049864, 0.012276154531876277,
    0.011833826239848241 };

  double zth_re;
  double zth_im;
  ac = fnu * 2.2250738585072014E-305;
  if ((fabs(z.re) <= ac) && (fabs(z.im) <= ac)) {
    zeta1->re = 1402.9773265065639 + fnu;
    zeta1->im = 0.0;
    zeta2->re = fnu;
    zeta2->im = 0.0;
    phi->re = 1.0;
    phi->im = 0.0;
    arg->re = 1.0;
    arg->im = 0.0;
  } else {
    zb.re = 1.0 / fnu * z.re;
    zb.im = 1.0 / fnu * z.im;
    ac = rt_powd_snf(fnu, 0.33333333333333331);
    fn23 = ac * ac;
    rfn13_re = 1.0 / ac;
    w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
    w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
    ac = rt_hypotd_snf(w2.re, w2.im);
    if (ac > 0.25) {
      b_sqrt(&w2);
      tmpr = w2.re;
      tmpi = w2.im;
      if (w2.re < 0.0) {
        tmpr = 0.0;
      }

      if (w2.im < 0.0) {
        tmpi = 0.0;
      }

      w2.re = tmpr;
      w2.im = tmpi;
      if (zb.im == 0.0) {
        if (tmpi == 0.0) {
          zc.re = (1.0 + tmpr) / zb.re;
          zc.im = 0.0;
        } else if (1.0 + tmpr == 0.0) {
          zc.re = 0.0;
          zc.im = tmpi / zb.re;
        } else {
          zc.re = (1.0 + tmpr) / zb.re;
          zc.im = tmpi / zb.re;
        }
      } else if (zb.re == 0.0) {
        if (1.0 + tmpr == 0.0) {
          zc.re = tmpi / zb.im;
          zc.im = 0.0;
        } else if (tmpi == 0.0) {
          zc.re = 0.0;
          zc.im = -((1.0 + tmpr) / zb.im);
        } else {
          zc.re = tmpi / zb.im;
          zc.im = -((1.0 + tmpr) / zb.im);
        }
      } else {
        brm = fabs(zb.re);
        ac = fabs(zb.im);
        if (brm > ac) {
          ac = zb.im / zb.re;
          pp = zb.re + ac * zb.im;
          zc.re = ((1.0 + tmpr) + ac * tmpi) / pp;
          zc.im = (tmpi - ac * (1.0 + tmpr)) / pp;
        } else if (ac == brm) {
          if (zb.re > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          if (zb.im > 0.0) {
            pp = 0.5;
          } else {
            pp = -0.5;
          }

          zc.re = ((1.0 + tmpr) * ac + tmpi * pp) / brm;
          zc.im = (tmpi * ac - (1.0 + tmpr) * pp) / brm;
        } else {
          ac = zb.re / zb.im;
          pp = zb.im + ac * zb.re;
          zc.re = (ac * (1.0 + tmpr) + tmpi) / pp;
          zc.im = (ac * tmpi - (1.0 + tmpr)) / pp;
        }
      }

      b_log(&zc);
      tmpr = zc.re;
      tmpi = zc.im;
      if (zc.re < 0.0) {
        tmpr = 0.0;
      }

      if (zc.im < 0.0) {
        tmpi = 0.0;
      } else {
        if (zc.im > 1.5707963267948966) {
          tmpi = 1.5707963267948966;
        }
      }

      zth_re = 1.5 * (tmpr - w2.re);
      zth_im = 1.5 * (tmpi - w2.im);
      zeta1->re = fnu * tmpr;
      zeta1->im = fnu * tmpi;
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      if ((zth_re >= 0.0) && (zth_im < 0.0)) {
        ac = 4.71238898038469;
      } else if (zth_re != 0.0) {
        ac = atan(zth_im / zth_re);
        if (zth_re < 0.0) {
          ac += 3.1415926535897931;
        }
      } else {
        ac = 1.5707963267948966;
      }

      pp = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
      ac *= 0.66666666666666663;
      zb.re = pp * cos(ac);
      zb.im = pp * sin(ac);
      if (zb.im < 0.0) {
        zb.im = 0.0;
      }

      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      if (zb.im == 0.0) {
        if (zth_im == 0.0) {
          tmpr = zth_re / zb.re;
          zth_im = 0.0;
        } else if (zth_re == 0.0) {
          tmpr = 0.0;
          zth_im /= zb.re;
        } else {
          tmpr = zth_re / zb.re;
          zth_im /= zb.re;
        }
      } else if (zb.re == 0.0) {
        if (zth_re == 0.0) {
          tmpr = zth_im / zb.im;
          zth_im = 0.0;
        } else if (zth_im == 0.0) {
          tmpr = 0.0;
          zth_im = -(zth_re / zb.im);
        } else {
          tmpr = zth_im / zb.im;
          zth_im = -(zth_re / zb.im);
        }
      } else {
        brm = fabs(zb.re);
        ac = zb.im;
        if (brm > ac) {
          ac = zb.im / zb.re;
          pp = zb.re + ac * zb.im;
          tmpr = (zth_re + ac * zth_im) / pp;
          zth_im = (zth_im - ac * zth_re) / pp;
        } else if (ac == brm) {
          if (zb.re > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          if (zb.im > 0.0) {
            pp = 0.5;
          } else {
            pp = -0.5;
          }

          tmpr = (zth_re * ac + zth_im * pp) / brm;
          zth_im = (zth_im * ac - zth_re * pp) / brm;
        } else {
          ac = zb.re / zb.im;
          pp = zb.im + ac * zb.re;
          tmpr = (ac * zth_re + zth_im) / pp;
          zth_im = (ac * zth_im - zth_re) / pp;
        }
      }

      if (w2.im == 0.0) {
        if (zth_im == 0.0) {
          suma.re = tmpr / w2.re;
          suma.im = 0.0;
        } else if (tmpr == 0.0) {
          suma.re = 0.0;
          suma.im = zth_im / w2.re;
        } else {
          suma.re = tmpr / w2.re;
          suma.im = zth_im / w2.re;
        }
      } else if (w2.re == 0.0) {
        if (tmpr == 0.0) {
          suma.re = zth_im / w2.im;
          suma.im = 0.0;
        } else if (zth_im == 0.0) {
          suma.re = 0.0;
          suma.im = -(tmpr / w2.im);
        } else {
          suma.re = zth_im / w2.im;
          suma.im = -(tmpr / w2.im);
        }
      } else {
        brm = fabs(w2.re);
        ac = fabs(w2.im);
        if (brm > ac) {
          ac = w2.im / w2.re;
          pp = w2.re + ac * w2.im;
          suma.re = (tmpr + ac * zth_im) / pp;
          suma.im = (zth_im - ac * tmpr) / pp;
        } else if (ac == brm) {
          if (w2.re > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          if (w2.im > 0.0) {
            pp = 0.5;
          } else {
            pp = -0.5;
          }

          suma.re = (tmpr * ac + zth_im * pp) / brm;
          suma.im = (zth_im * ac - tmpr * pp) / brm;
        } else {
          ac = w2.re / w2.im;
          pp = w2.im + ac * w2.re;
          suma.re = (ac * tmpr + zth_im) / pp;
          suma.im = (ac * zth_im - tmpr) / pp;
        }
      }

      zb.re = suma.re + suma.re;
      zb.im = suma.im + suma.im;
      b_sqrt(&zb);
      phi->re = zb.re * rfn13_re - zb.im * 0.0;
      phi->im = zb.re * 0.0 + zb.im * rfn13_re;
    } else {
      memset(&p[0], 0, 30U * sizeof(creal_T));
      memset(&ap[0], 0, 30U * sizeof(double));
      p[0].re = 1.0;
      p[0].im = 0.0;
      suma.re = 0.6299605249474366;
      suma.im = 0.0;
      ap[0] = 1.0;
      if (ac >= tol) {
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i + 1 < 31)) {
          pp = p[i - 1].re;
          p[i].re = p[i - 1].re * w2.re - p[i - 1].im * w2.im;
          p[i].im = pp * w2.im + p[i - 1].im * w2.re;
          suma.re += p[i].re * GAMA[i];
          suma.im += p[i].im * GAMA[i];
          ap[i] = ap[i - 1] * ac;
          if (ap[i] < tol) {
            exitg1 = true;
          } else {
            i++;
          }
        }
      }

      zb.re = w2.re * suma.re - w2.im * suma.im;
      zb.im = w2.re * suma.im + w2.im * suma.re;
      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      b_sqrt(&suma);
      b_sqrt(&w2);
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      ac = (zb.re * suma.re - zb.im * suma.im) * 0.66666666666666663 + 1.0;
      pp = (zb.re * suma.im + zb.im * suma.re) * 0.66666666666666663;
      zeta1->re = zeta2->re * ac - zeta2->im * pp;
      zeta1->im = zeta2->re * pp + zeta2->im * ac;
      suma.re += suma.re;
      suma.im += suma.im;
      b_sqrt(&suma);
      phi->re = suma.re * rfn13_re - suma.im * 0.0;
      phi->im = suma.re * 0.0 + suma.im * rfn13_re;
    }
  }
}

/*
 * File trailer for cunhj.c
 *
 * [EOF]
 */
