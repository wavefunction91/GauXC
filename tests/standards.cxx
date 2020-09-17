#include "standards.hpp"
#include "ut_common.hpp"
#include "basis/parse_basis.hpp"

namespace GauXC {

Molecule make_water() {

  Molecule mol;
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028,  2.174611055780858);
  mol.emplace_back(AtomicNumber(8), 0., 0.000000000000000,  0.000000000000000);
  mol.emplace_back(AtomicNumber(1), 0., 1.579252144093028, -2.174611055780858);

  return mol;

}

Molecule make_benzene() {

  Molecule mol;
  mol.emplace_back(AtomicNumber(6),  6.92768e-01,  -1.77656e+00,   1.40218e-03);
  mol.emplace_back(AtomicNumber(6),  3.35108e+00,  -1.77668e+00,   2.21098e-03);
  mol.emplace_back(AtomicNumber(6),  4.68035e+00,   5.25219e-01,   1.22454e-03);
  mol.emplace_back(AtomicNumber(6),  3.35121e+00,   2.82744e+00,  -7.02978e-04);
  mol.emplace_back(AtomicNumber(6),  6.93087e-01,   2.82756e+00,  -1.55902e-03);
  mol.emplace_back(AtomicNumber(6), -6.36278e-01,   5.25491e-01,  -4.68652e-04);
  mol.emplace_back(AtomicNumber(1), -3.41271e-01,  -3.56759e+00,   2.21287e-03);
  mol.emplace_back(AtomicNumber(1),  4.38492e+00,  -3.56783e+00,   3.73599e-03);
  mol.emplace_back(AtomicNumber(1),  6.74844e+00,   5.25274e-01,   1.88028e-03);
  mol.emplace_back(AtomicNumber(1),  4.38551e+00,   4.61832e+00,  -1.48721e-03);
  mol.emplace_back(AtomicNumber(1), -3.41001e-01,   4.61857e+00,  -3.05569e-03);
  mol.emplace_back(AtomicNumber(1), -2.70437e+00,   5.25727e-01,  -1.09793e-03);

  return mol;
}

Molecule make_taxol() {

  Molecule mol;
#if 0
  mol.emplace_back(AtomicNumber(6), -6.3960615647,   1.6375500661,  -1.0167390582);
  mol.emplace_back(AtomicNumber(6), -5.5976011922,   1.9479656072,    .3156763742);
  mol.emplace_back(AtomicNumber(6), -4.1358423840,   1.3888885892,    .3061559180);
  mol.emplace_back(AtomicNumber(6), -4.0900116402,   -.2069339772,    .0759908147);
  mol.emplace_back(AtomicNumber(6), -5.3005244050,   -.5559695927,   -.7831189443);
  mol.emplace_back(AtomicNumber(6), -5.6408352846,    .5560151731,  -1.7454174478);
  mol.emplace_back(AtomicNumber(6), -2.8203119095,   -.7723427440,   -.5951161079);
  mol.emplace_back(AtomicNumber(6), -3.4405167500,   1.9809801074,   1.5954932385);
  mol.emplace_back(AtomicNumber(6), -1.5722954629,   -.0769638917,   -.1719253097);
  mol.emplace_back(AtomicNumber(8), -2.6953275274,  -2.1590273791,   -.1697708437);
  mol.emplace_back(AtomicNumber(8), -1.8693493295,   1.2308626735,   -.7751223509);
  mol.emplace_back(AtomicNumber(6), -1.3334064307,   -.1176307913,   1.4101576262);
  mol.emplace_back(AtomicNumber(6),  -.8881044770,   1.2874584739,   1.4845316346);
  mol.emplace_back(AtomicNumber(6),   .3196036852,   1.5770567320,    .9640605315);
  mol.emplace_back(AtomicNumber(8), -4.0360627508,   2.4497903759,   2.5047989338);
  mol.emplace_back(AtomicNumber(8), -1.9504980519,   2.1870799707,   1.7492357544);
  mol.emplace_back(AtomicNumber(6),  -.2198636798,   -.3744036919,   -.8072269875);
  mol.emplace_back(AtomicNumber(6),   .8361321957,    .5392960246,   -.0882305557);
  mol.emplace_back(AtomicNumber(8), -6.3705650302,   -.7461436866,    .1814962812);
  mol.emplace_back(AtomicNumber(6), -7.5079885256,  -1.5312544379,    .1764442320);
  mol.emplace_back(AtomicNumber(8), -7.7343168261,  -2.1872429459,   1.1532459677);
  mol.emplace_back(AtomicNumber(6), -8.4523249942,  -1.4157500320,   -.9939259615);
  mol.emplace_back(AtomicNumber(6), -5.4161056771,  -1.5091023683,  -2.0070002809);
  mol.emplace_back(AtomicNumber(8), -6.2192571555,   -.4454467295,  -2.6902950026);
  mol.emplace_back(AtomicNumber(8), -5.4267626184,   3.3774437107,    .3776243498);
  mol.emplace_back(AtomicNumber(8), -1.6122441118,   4.1801170473,    .0343779432);
  mol.emplace_back(AtomicNumber(6),  -.7940923698,   5.2257739597,    .2990567126);
  mol.emplace_back(AtomicNumber(8),  -.3048947277,   5.8645765943,   -.5895332555);
  mol.emplace_back(AtomicNumber(6),  -.5432141747,   5.4736984394,   1.7763873918);
  mol.emplace_back(AtomicNumber(6),  1.2298580068,   2.7852908502,   1.0315034732);
  mol.emplace_back(AtomicNumber(8),  1.9361230485,   -.3334212217,    .3599506557);
  mol.emplace_back(AtomicNumber(6),  3.1923971442,    .0392867864,    .0732007951);
  mol.emplace_back(AtomicNumber(8),  3.4564612550,   1.0474089559,   -.5229789642);
  mol.emplace_back(AtomicNumber(6),  4.1741649865,  -1.0074733119,    .5409360653);
  mol.emplace_back(AtomicNumber(6),  5.5923102950,   -.5507792057,    .2305424754);
  mol.emplace_back(AtomicNumber(6),  6.6535032780,  -1.6276654266,    .4267849964);
  mol.emplace_back(AtomicNumber(6),  6.3531145194,  -2.9590427132,    .6578478411);
  mol.emplace_back(AtomicNumber(6),  7.3484376981,  -3.9223474323,    .5946432758);
  mol.emplace_back(AtomicNumber(6),  8.6553264573,  -3.5577622028,    .3262773335);
  mol.emplace_back(AtomicNumber(6),  8.9721236891,  -2.2174054241,    .1751067514);
  mol.emplace_back(AtomicNumber(6),  7.9780880254,  -1.2615408527,    .2323472017);
  mol.emplace_back(AtomicNumber(8),  4.0117554647,  -1.2781153298,   1.9325731679);
  mol.emplace_back(AtomicNumber(7),  5.9708266751,    .6810003374,    .8873415746);
  mol.emplace_back(AtomicNumber(6),  6.1083010659,   1.9571429450,    .3956141999);
  mol.emplace_back(AtomicNumber(8),  6.3193157933,   2.8803620623,   1.1682891449);
  mol.emplace_back(AtomicNumber(6),  5.8856236342,   2.2294771536,  -1.0684368500);
  mol.emplace_back(AtomicNumber(6),  6.1182756121,   1.3530945943,  -2.1154261935);
  mol.emplace_back(AtomicNumber(6),  5.5666092653,   1.5985926913,  -3.3627521786);
  mol.emplace_back(AtomicNumber(6),  4.8504992902,   2.7660991077,  -3.5880645167);
  mol.emplace_back(AtomicNumber(6),  4.7600523715,   3.7102511809,  -2.5859738609);
  mol.emplace_back(AtomicNumber(6),  5.2758308332,   3.4419993713,  -1.3326644938);
  mol.emplace_back(AtomicNumber(6), -2.1366949616,  -3.0476305872,  -1.0300505942);
  mol.emplace_back(AtomicNumber(8), -2.1451655154,  -2.8533797614,  -2.2328870532);
  mol.emplace_back(AtomicNumber(6), -1.5033714554,  -4.2030087508,   -.3691075134);
  mol.emplace_back(AtomicNumber(6), -1.4395726939,  -4.3266981108,   1.0054185593);
  mol.emplace_back(AtomicNumber(6),  -.7638829811,  -5.3989751504,   1.5627241997);
  mol.emplace_back(AtomicNumber(6),  -.1728504183,  -6.3468606504,    .7427716273);
  mol.emplace_back(AtomicNumber(6),  -.2555829630,  -6.2275386962,   -.6353203163);
  mol.emplace_back(AtomicNumber(6),  -.9155593963,  -5.1509643641,  -1.1911103638);
  mol.emplace_back(AtomicNumber(6),  -.2097392796,  -1.1547474452,   1.8364340368);
  mol.emplace_back(AtomicNumber(6), -2.4300502071,   -.4562312768,   2.4707641203);
  mol.emplace_back(AtomicNumber(1), -3.6491742439,   1.8683152472,   -.5155830220);
  mol.emplace_back(AtomicNumber(1), -4.2442039674,   -.7298028968,    .9946070928);
  mol.emplace_back(AtomicNumber(1), -6.1319126304,   1.5761651125,   1.1739416597);
  mol.emplace_back(AtomicNumber(1), -6.4278095287,   2.5627462321,  -1.5682088644);
  mol.emplace_back(AtomicNumber(1), -7.3896383806,   1.3053381318,   -.7744956635);
  mol.emplace_back(AtomicNumber(1), -4.7465705868,    .9454137315,  -2.2120480815);
  mol.emplace_back(AtomicNumber(1), -2.8855554103,   -.7190312325,  -1.6682581370);
  mol.emplace_back(AtomicNumber(1), -1.1756471657,   1.8958446516,   -.6850815225);
  mol.emplace_back(AtomicNumber(1), -1.8481632425,   3.6442035335,    .8091515684);
  mol.emplace_back(AtomicNumber(1),  -.3040222781,   -.0416779115,  -1.8339005061);
  mol.emplace_back(AtomicNumber(1),   .0915191278,  -1.3952921601,   -.8028715623);
  mol.emplace_back(AtomicNumber(1),  1.2665207958,   1.1759763221,   -.8367612418);
  mol.emplace_back(AtomicNumber(1), -7.9447751635,  -1.4087835985,  -1.9396806789);
  mol.emplace_back(AtomicNumber(1), -8.9860724820,   -.4736978340,   -.9254124624);
  mol.emplace_back(AtomicNumber(1), -9.1714803590,  -2.2205724001,   -.9250166491);
  mol.emplace_back(AtomicNumber(1), -5.9696659013,  -2.4206608341,  -1.8715269767);
  mol.emplace_back(AtomicNumber(1), -4.4956987862,  -1.6916167197,  -2.5341522638);
  mol.emplace_back(AtomicNumber(1), -6.1859306778,   3.7990292112,    .8057985278);
  mol.emplace_back(AtomicNumber(1),  -.4047348600,   4.5302819326,   2.2918679722);
  mol.emplace_back(AtomicNumber(1), -1.4092253680,   5.9590265117,   2.2201988007);
  mol.emplace_back(AtomicNumber(1),   .3171040869,   6.0931663404,   1.8871608880);
  mol.emplace_back(AtomicNumber(1),  1.4832824150,   3.0343295148,   2.0569031127);
  mol.emplace_back(AtomicNumber(1),   .8689347198,   3.6687756404,    .5365846962);
  mol.emplace_back(AtomicNumber(1),  2.1684888568,   2.5341563953,    .5367744746);
  mol.emplace_back(AtomicNumber(1),  3.9669613366,  -1.8924968836,   -.0572581027);
  mol.emplace_back(AtomicNumber(1),  5.5541772134,   -.3516751663,   -.8217997298);
  mol.emplace_back(AtomicNumber(1),  5.3544974548,  -3.2557099687,    .9038847941);
  mol.emplace_back(AtomicNumber(1),  7.0989139525,  -4.9513958159,    .7469287918);
  mol.emplace_back(AtomicNumber(1),  9.4171939156,  -4.3036934271,    .2407856159);
  mol.emplace_back(AtomicNumber(1),  9.9857443995,  -1.9177653222,    .0023467668);
  mol.emplace_back(AtomicNumber(1),  8.2254733522,   -.2253064772,    .1253413078);
  mol.emplace_back(AtomicNumber(1),  3.1299494847,  -1.6337868182,   2.1083778803);
  mol.emplace_back(AtomicNumber(1),  6.0936138326,    .6250743708,   1.8859547914);
  mol.emplace_back(AtomicNumber(1),  6.7084379892,    .4705180459,  -1.9833352247);
  mol.emplace_back(AtomicNumber(1),  5.6743950704,    .8727859493,  -4.1412749888);
  mol.emplace_back(AtomicNumber(1),  4.3798456238,   2.9274160428,  -4.5335955847);
  mol.emplace_back(AtomicNumber(1),  4.2563843627,   4.6442138650,  -2.7660687217);
  mol.emplace_back(AtomicNumber(1),  5.1747305805,   4.1495326330,   -.5327311357);
  mol.emplace_back(AtomicNumber(1), -1.9015964720,  -3.5961881882,   1.6316249866);
  mol.emplace_back(AtomicNumber(1),  -.6959960683,  -5.4955297144,   2.6284271470);
  mol.emplace_back(AtomicNumber(1),   .3482320371,  -7.1761491687,   1.1766068045);
  mol.emplace_back(AtomicNumber(1),   .1934856075,  -6.9647950021,  -1.2689411336);
  mol.emplace_back(AtomicNumber(1),  -.9826070096,  -5.0258752831,  -2.2518209168);
  mol.emplace_back(AtomicNumber(1),   .6895305096,   -.6503325497,   2.1318038533);
  mol.emplace_back(AtomicNumber(1),   .0408559146,  -1.8507203279,   1.0570933960);
  mol.emplace_back(AtomicNumber(1),  -.5497152190,  -1.7387966223,   2.6730562622);
  mol.emplace_back(AtomicNumber(1), -2.6864693812,  -1.5071204618,   2.3903167090);
  mol.emplace_back(AtomicNumber(1), -3.3390142876,    .1000510598,   2.4688669158);
  mol.emplace_back(AtomicNumber(1), -1.9770019524,   -.2959025420,   3.4410100992);
  
  // Convert to Bohr
  for( auto& atom : mol ) {
    atom.x *= 1.8897161646321;
    atom.y *= 1.8897161646321;
    atom.z *= 1.8897161646321;
  }
#else
 mol.emplace_back(AtomicNumber(6), -1.2086803766057596e+01,  3.0945209181978379e+00, -1.9213582223119234e+00);
 mol.emplace_back(AtomicNumber(6), -1.0577932448957723e+01,  3.6811212336040056e+00,  5.9654184843902913e-01);
 mol.emplace_back(AtomicNumber(6), -7.8156088394525174e+00,  2.6246188628367846e+00,  5.7855079493075290e-01);
 mol.emplace_back(AtomicNumber(6), -7.7290012917984576e+00, -3.9104851472197344e-01,  1.4360181746387324e-01);
 mol.emplace_back(AtomicNumber(6), -1.0016538723457261e+01, -1.0506301884189346e+00, -1.4798802215219533e+00);
 mol.emplace_back(AtomicNumber(6), -1.0659633036976830e+01,  1.0507163228854037e+00, -3.2983607127617107e+00);
 mol.emplace_back(AtomicNumber(6), -5.3296167124683658e+00, -1.4595161557523739e+00, -1.1246063755711582e+00);
 mol.emplace_back(AtomicNumber(6), -6.5016339180648162e+00,  3.7435095926457911e+00,  3.0150450380672256e+00);
 mol.emplace_back(AtomicNumber(6), -2.9712075986289155e+00, -1.4544066636007139e-01, -3.2489172590696380e-01);
 mol.emplace_back(AtomicNumber(8), -5.0934304773948895e+00, -4.0799701492478251e+00, -3.2082037551434694e-01);
 mol.emplace_back(AtomicNumber(8), -3.5325580104758747e+00,  2.3259931830029714e+00, -1.4647688511505077e+00);
 mol.emplace_back(AtomicNumber(6), -2.5197727859935175e+00, -2.2228996342624507e-01,  2.6648115148166873e+00);
 mol.emplace_back(AtomicNumber(6), -1.6782741111341526e+00,  2.4329437378871082e+00,  2.8053580113962715e+00);
 mol.emplace_back(AtomicNumber(6),  6.0396339010261468e-01,  2.9802050925878079e+00,  1.8218102413447033e+00);
 mol.emplace_back(AtomicNumber(8), -7.6270526734215913e+00,  4.6294325409403099e+00,  4.7333836424213507e+00);
 mol.emplace_back(AtomicNumber(8), -3.6859068601693008e+00,  4.1329818606531488e+00,  3.3055762659777010e+00);
 mol.emplace_back(AtomicNumber(6), -4.1548210975523431e-01, -7.0752038696097874e-01, -1.5254378173009282e+00);
 mol.emplace_back(AtomicNumber(6),  1.5800607404539240e+00,  1.0191217134510033e+00, -1.6673157413020209e-01);
 mol.emplace_back(AtomicNumber(8), -1.2038622302183510e+01, -1.4100071160962910e+00,  3.4297823949049211e-01);
 mol.emplace_back(AtomicNumber(6), -1.4188041041940112e+01, -2.8936513070712166e+00,  3.3343125081954544e-01);
 mol.emplace_back(AtomicNumber(8), -1.4615739513441163e+01, -4.1332898391241510e+00,  2.1793188768721445e+00);
 mol.emplace_back(AtomicNumber(6), -1.5972578209014014e+01, -2.6753796293979817e+00, -1.8782477205883634e+00);
 mol.emplace_back(AtomicNumber(6), -1.0234955657186312e+01, -2.8517899654379599e+00, -3.7926805907470298e+00);
 mol.emplace_back(AtomicNumber(8), -1.1752691879022565e+01, -8.4177226145120299e-01, -5.0839203844900425e+00);
 mol.emplace_back(AtomicNumber(8), -1.0255094356124170e+01,  6.3824431564943875e+00,  7.1360654789628697e-01);
 mol.emplace_back(AtomicNumber(8), -3.0466995986806813e+00,  7.8992758213447534e+00,  6.4964892713405817e-02);
 mol.emplace_back(AtomicNumber(6), -1.5006169888776586e+00,  9.8752808642845284e+00,  5.6513524198512377e-01);
 mol.emplace_back(AtomicNumber(8), -5.7616749084376817e-01,  1.1082442804729819e+01, -1.1140563142981272e+00);
 mol.emplace_back(AtomicNumber(6), -1.0265259435237764e+00,  1.0343790196882921e+01,  3.3568854208163854e+00);
 mol.emplace_back(AtomicNumber(6),  2.3240946382296985e+00,  5.2634365065468458e+00,  1.9492589210498048e+00);
 mol.emplace_back(AtomicNumber(8),  3.6587420426523574e+00, -6.3007474793062079e-01,  6.8020810883388105e-01);
 mol.emplace_back(AtomicNumber(6),  6.0327558506041203e+00,  7.4241261284371751e-02,  1.3832944491593385e-01);
 mol.emplace_back(AtomicNumber(8),  6.5317646635450561e+00,  1.9793159250755847e+00, -9.8828694034904063e-01);
 mol.emplace_back(AtomicNumber(6),  7.8880280573628845e+00, -1.9038485007213328e+00,  1.0222209409848111e+00);
 mol.emplace_back(AtomicNumber(6),  1.0567934103013757e+01, -1.0408217792120671e+00,  4.3566210733177319e-01);
 mol.emplace_back(AtomicNumber(6),  1.2573298062333292e+01, -3.0758416580427919e+00,  8.0650669941235142e-01);
 mol.emplace_back(AtomicNumber(6),  1.2005645618403424e+01, -5.5917799176951135e+00,  1.2431521621342125e+00);
 mol.emplace_back(AtomicNumber(6),  1.3886533696646906e+01, -7.4121618807047280e+00,  1.1237128524633548e+00);
 mol.emplace_back(AtomicNumber(6),  1.6356195349639108e+01, -6.7231956973130487e+00,  6.1657475673657036e-01);
 mol.emplace_back(AtomicNumber(6),  1.6954855311814825e+01, -4.1902886580713368e+00,  3.3090377896994216e-01);
 mol.emplace_back(AtomicNumber(6),  1.5076400284128072e+01, -2.3839665355324109e+00,  4.3907254552391500e-01);
 mol.emplace_back(AtomicNumber(8),  7.5811185631563625e+00, -2.4152877556623662e+00,  3.6520337410246904e+00);
 mol.emplace_back(AtomicNumber(7),  1.1283226343750929e+01,  1.2869040361025486e+00,  1.6768324346418022e+00);
 mol.emplace_back(AtomicNumber(6),  1.1543015272867631e+01,  3.6984638873544977e+00,  7.4760243516847125e-01);
 mol.emplace_back(AtomicNumber(8),  1.1941775287297162e+01,  5.4430950468579473e+00,  2.2077463597841165e+00);
 mol.emplace_back(AtomicNumber(6),  1.1122215943020368e+01,  4.2131009190396647e+00, -2.0190528830502950e+00);
 mol.emplace_back(AtomicNumber(6),  1.1561864432050251e+01,  2.5569780204241210e+00, -3.9975758556682925e+00);
 mol.emplace_back(AtomicNumber(6),  1.0519366199245606e+01,  3.0209021545750638e+00, -6.3546801864667897e+00);
 mol.emplace_back(AtomicNumber(6),  9.1661145683169938e+00,  5.2271693719704002e+00, -6.7804587674167145e+00);
 mol.emplace_back(AtomicNumber(6),  8.9951946754246315e+00,  7.0113580822646702e+00, -4.8867820118174006e+00);
 mol.emplace_back(AtomicNumber(6),  9.9698746390655639e+00,  6.5044356660672706e+00, -2.5183707285513894e+00);
 mol.emplace_back(AtomicNumber(6), -4.0377679995008773e+00, -5.7591867255031710e+00, -1.9465133778446326e+00);
 mol.emplace_back(AtomicNumber(8), -4.0537750251579601e+00, -5.3921058916041984e+00, -4.2195446949336652e+00);
 mol.emplace_back(AtomicNumber(6), -2.8409601103901343e+00, -7.9425348683811841e+00, -6.9751206080714567e-01);
 mol.emplace_back(AtomicNumber(6), -2.7203979327175718e+00, -8.1762738665359613e+00,  1.8999655813321477e+00);
 mol.emplace_back(AtomicNumber(6), -1.4435295219394657e+00, -1.0202583655676063e+01,  2.9531205338123163e+00);
 mol.emplace_back(AtomicNumber(6), -3.2663992767103123e-01, -1.1993827519622323e+01,  1.4036348480006320e+00);
 mol.emplace_back(AtomicNumber(6), -4.8298176752672539e-01, -1.1768341721712316e+01, -1.2005813130518102e+00);
 mol.emplace_back(AtomicNumber(6), -1.7301563856612605e+00, -9.7339112272526283e+00, -2.2508722102401046e+00);
 mol.emplace_back(AtomicNumber(6), -3.9634976757425749e-01, -2.1821562579257932e+00,  3.4703571264251427e+00);
 mol.emplace_back(AtomicNumber(6), -4.5921290309317024e+00, -8.6215210076361271e-01,  4.6690671708196323e+00);
 mol.emplace_back(AtomicNumber(1), -6.8959394070872548e+00,  3.5306038782787996e+00, -9.7431063616055869e-01);
 mol.emplace_back(AtomicNumber(1), -8.0203825398126884e+00, -1.3791275009304449e+00,  1.8795348721078948e+00);
 mol.emplace_back(AtomicNumber(1), -1.1587634659944230e+01,  2.9785201760463589e+00,  2.2184280639048839e+00);
 mol.emplace_back(AtomicNumber(1), -1.2146798718726231e+01,  4.8428881580111955e+00, -2.9634850472368570e+00);
 mol.emplace_back(AtomicNumber(1), -1.3964391697131694e+01,  2.4667313920951672e+00, -1.4635845836837487e+00);
 mol.emplace_back(AtomicNumber(1), -8.9697177964989407e+00,  1.7865728987730178e+00, -4.1801647485281395e+00);
 mol.emplace_back(AtomicNumber(1), -5.4529090515434682e+00, -1.3587720069579514e+00, -3.1525507578496224e+00);
 mol.emplace_back(AtomicNumber(1), -2.2216510029174796e+00,  3.5826269092351706e+00, -1.2946163576519383e+00);
 mol.emplace_back(AtomicNumber(1), -3.4925221112667595e+00,  6.8865461264605816e+00,  1.5290747478455913e+00);
 mol.emplace_back(AtomicNumber(1), -5.7451880016055556e-01, -7.8759832528791970e-02, -3.4655694476174230e+00);
 mol.emplace_back(AtomicNumber(1),  1.7294607429427239e-01, -2.6367198571889188e+00, -1.5172072571073427e+00);
 mol.emplace_back(AtomicNumber(1),  2.3933772634322219e+00,  2.2222730183210051e+00, -1.5812494652173732e+00);
 mol.emplace_back(AtomicNumber(1), -1.5013448103227674e+01, -2.6622149789623912e+00, -3.6654649892784938e+00);
 mol.emplace_back(AtomicNumber(1), -1.6981214708273136e+01, -8.9515910784280783e-01, -1.7487759807417653e+00);
 mol.emplace_back(AtomicNumber(1), -1.7331584792005351e+01, -4.1962733749250765e+00, -1.7480280020619636e+00);
 mol.emplace_back(AtomicNumber(1), -1.1281032799333719e+01, -4.5743856887531873e+00, -3.5366731669845874e+00);
 mol.emplace_back(AtomicNumber(1), -8.4956388349978944e+00, -3.1966920786440185e+00, -4.7888533929860433e+00);
 mol.emplace_back(AtomicNumber(1), -1.1689713967991045e+01,  7.1791242333748100e+00,  1.5227384198815992e+00);
 mol.emplace_back(AtomicNumber(1), -7.6483798359627653e-01,  8.5609915055313675e+00,  4.3310024704230692e+00);
 mol.emplace_back(AtomicNumber(1), -2.6630498022676887e+00,  1.1260927268299502e+01,  4.1955673744294213e+00);
 mol.emplace_back(AtomicNumber(1),  5.9923983423304450e-01,  1.1514414788753900e+01,  3.5662169754779183e+00);
 mol.emplace_back(AtomicNumber(1),  2.8029973286521832e+00,  5.7340513433073204e+00,  3.8869832689241863e+00);
 mol.emplace_back(AtomicNumber(1),  1.6420485227504928e+00,  6.9329806754739982e+00,  1.0139980457088096e+00);
 mol.emplace_back(AtomicNumber(1),  4.0978497495518589e+00,  4.7888612003889675e+00,  1.0143566748834405e+00);
 mol.emplace_back(AtomicNumber(1),  7.4964699351311967e+00, -3.5763005450404277e+00, -1.0820212475302107e-01);
 mol.emplace_back(AtomicNumber(1),  1.0495873027673579e+01, -6.6456970144300698e-01, -1.5529763071562379e+00);
 mol.emplace_back(AtomicNumber(1),  1.0118532998369913e+01, -6.1523997404987671e+00,  1.7080945864726838e+00);
 mol.emplace_back(AtomicNumber(1),  1.3415002189713961e+01, -9.3567813551320889e+00,  1.4114907497968301e+00);
 mol.emplace_back(AtomicNumber(1),  1.7795916085761995e+01, -8.1328013178793466e+00,  4.5501883614360161e-01);
 mol.emplace_back(AtomicNumber(1),  1.8870320711246347e+01, -3.6240509701642987e+00,  4.4347462120823654e-03);
 mol.emplace_back(AtomicNumber(1),  1.5543890765479290e+01, -4.2576750545487596e-01,  2.3686072684490841e-01);
 mol.emplace_back(AtomicNumber(1),  5.9147468854947478e+00, -3.0874094108381582e+00,  3.9842564750356408e+00);
 mol.emplace_back(AtomicNumber(1),  1.1515260426394116e+01,  1.1812192835585826e+00,  3.5639377833876540e+00);
 mol.emplace_back(AtomicNumber(1),  1.2677109613786142e+01,  8.8915017963072496e-01, -3.7479601190147447e+00);
 mol.emplace_back(AtomicNumber(1),  1.0723051836388365e+01,  1.6493262912262463e+00, -7.8258749739310440e+00);
 mol.emplace_back(AtomicNumber(1),  8.2767081031027754e+00,  5.5320141766946964e+00, -8.5672534000232403e+00);
 mol.emplace_back(AtomicNumber(1),  8.0434001493673932e+00,  8.7762916391646382e+00, -5.2271119507564983e+00);
 mol.emplace_back(AtomicNumber(1),  9.7788228640439065e+00,  7.8414796587836983e+00, -1.0067158722817757e+00);
 mol.emplace_back(AtomicNumber(1), -3.5934962737291110e+00, -6.7958102805763634e+00,  3.0833241414797969e+00);
 mol.emplace_back(AtomicNumber(1), -1.3152418585083292e+00, -1.0385045324623427e+01,  4.9670070898790231e+00);
 mol.emplace_back(AtomicNumber(1),  6.5806313071028222e-01, -1.3560955585033136e+01,  2.2234644572978923e+00);
 mol.emplace_back(AtomicNumber(1),  3.6563478099020330e-01, -1.3161554123525679e+01, -2.3979510386750413e+00);
 mol.emplace_back(AtomicNumber(1), -1.8568580030146926e+00, -9.4975271399468024e+00, -4.2553245090507668e+00);
 mol.emplace_back(AtomicNumber(1),  1.3030237241995339e+00, -1.2289503206607242e+00,  4.0285251450313533e+00);
 mol.emplace_back(AtomicNumber(1),  7.7206483623984537e-02, -3.4973543020032318e+00,  1.9976168632214686e+00);
 mol.emplace_back(AtomicNumber(1), -1.0388111358931265e+00, -3.2858491667457268e+00,  5.0513438887385380e+00);
 mol.emplace_back(AtomicNumber(1), -5.0766910083063879e+00, -2.8480447052171418e+00,  4.5170436069382500e+00);
 mol.emplace_back(AtomicNumber(1), -6.3098220769200406e+00,  1.8906908793105312e-01,  4.6654819741695350e+00);
 mol.emplace_back(AtomicNumber(1), -3.7359919697540209e+00, -5.5917472382856404e-01,  6.5025662128697084e+00);

#endif
  return mol;
}

Molecule make_ubiquitin() {

  Molecule mol;

#if 0
  mol.emplace_back(AtomicNumber(7 ), 27.34000000,  24.43000000,  2.61400000);
  mol.emplace_back(AtomicNumber(6 ), 26.26600000,  25.41300000,  2.84200000);
  mol.emplace_back(AtomicNumber(6 ), 26.91300000,  26.63900000,  3.53100000);
  mol.emplace_back(AtomicNumber(8 ), 27.88600000,  26.46300000,  4.26300000);
  mol.emplace_back(AtomicNumber(6 ), 25.11200000,  24.88000000,  3.64900000);
  mol.emplace_back(AtomicNumber(6 ), 25.35300000,  24.86000000,  5.13400000);
  mol.emplace_back(AtomicNumber(16), 23.93000000,  23.95900000,  5.90400000);
  mol.emplace_back(AtomicNumber(6 ), 24.44700000,  23.98400000,  7.62000000);
  mol.emplace_back(AtomicNumber(7 ), 26.33500000,  27.77000000,  3.25800000);
  mol.emplace_back(AtomicNumber(6 ), 26.85000000,  29.02100000,  3.89800000);
  mol.emplace_back(AtomicNumber(6 ), 26.10000000,  29.25300000,  5.20200000);
  mol.emplace_back(AtomicNumber(8 ), 24.86500000,  29.02400000,  5.33000000);
  mol.emplace_back(AtomicNumber(6 ), 26.73300000,  30.14800000,  2.90500000);
  mol.emplace_back(AtomicNumber(6 ), 26.88200000,  31.54600000,  3.40900000);
  mol.emplace_back(AtomicNumber(6 ), 26.78600000,  32.56200000,  2.27000000);
  mol.emplace_back(AtomicNumber(8 ), 27.78300000,  33.16000000,  1.87000000);
  mol.emplace_back(AtomicNumber(7 ), 25.56200000,  32.73300000,  1.80600000);
  mol.emplace_back(AtomicNumber(7 ), 26.84900000,  29.65600000,  6.21700000);
  mol.emplace_back(AtomicNumber(6 ), 26.23500000,  30.05800000,  7.49700000);
  mol.emplace_back(AtomicNumber(6 ), 26.88200000,  31.42800000,  7.86200000);
  mol.emplace_back(AtomicNumber(8 ), 27.90600000,  31.71100000,  7.26400000);
  mol.emplace_back(AtomicNumber(6 ), 26.34400000,  29.05000000,  8.64500000);
  mol.emplace_back(AtomicNumber(6 ), 27.81000000,  28.74800000,  8.99900000);
  mol.emplace_back(AtomicNumber(6 ), 25.49100000,  27.77100000,  8.28700000);
  mol.emplace_back(AtomicNumber(6 ), 27.96700000,  28.08700000, 10.41700000);
  mol.emplace_back(AtomicNumber(7 ), 26.21400000,  32.09700000,  8.77100000);
  mol.emplace_back(AtomicNumber(6 ), 26.77200000,  33.43600000,  9.19700000);
  mol.emplace_back(AtomicNumber(6 ), 27.15100000,  33.36200000, 10.65000000);
  mol.emplace_back(AtomicNumber(8 ), 26.35000000,  32.77800000, 11.39500000);
  mol.emplace_back(AtomicNumber(6 ), 25.69500000,  34.49800000,  8.94600000);
  mol.emplace_back(AtomicNumber(6 ), 25.28800000,  34.60900000,  7.49900000);
  mol.emplace_back(AtomicNumber(6 ), 24.14700000,  33.96600000,  7.03800000);
  mol.emplace_back(AtomicNumber(6 ), 26.13600000,  35.34600000,  6.64000000);
  mol.emplace_back(AtomicNumber(6 ), 23.81200000,  34.03100000,  5.67700000);
  mol.emplace_back(AtomicNumber(6 ), 25.81000000,  35.39200000,  5.26700000);
  mol.emplace_back(AtomicNumber(6 ), 24.62000000,  34.77800000,  4.85300000);
  mol.emplace_back(AtomicNumber(7 ), 28.26000000,  33.94300000, 11.09600000);
  mol.emplace_back(AtomicNumber(6 ), 28.60500000,  33.96500000, 12.50300000);
  mol.emplace_back(AtomicNumber(6 ), 28.63800000,  35.46100000, 12.90000000);
  mol.emplace_back(AtomicNumber(8 ), 29.52200000,  36.10300000, 12.32000000);
  mol.emplace_back(AtomicNumber(6 ), 29.96300000,  33.31700000, 12.81400000);
  mol.emplace_back(AtomicNumber(6 ), 30.21100000,  33.39400000, 14.30400000);
  mol.emplace_back(AtomicNumber(6 ), 29.95700000,  31.83800000, 12.35200000);
  mol.emplace_back(AtomicNumber(7 ), 27.75100000,  35.86700000, 13.74000000);
  mol.emplace_back(AtomicNumber(6 ), 27.69100000,  37.31500000, 14.14300000);
  mol.emplace_back(AtomicNumber(6 ), 28.46900000,  37.47500000, 15.42000000);
  mol.emplace_back(AtomicNumber(8 ), 28.21300000,  36.75300000, 16.41100000);
  mol.emplace_back(AtomicNumber(6 ), 26.21900000,  37.68400000, 14.30700000);
  mol.emplace_back(AtomicNumber(6 ), 25.88400000,  39.13900000, 14.61500000);
  mol.emplace_back(AtomicNumber(6 ), 24.34800000,  39.29600000, 14.64200000);
  mol.emplace_back(AtomicNumber(6 ), 23.86500000,  40.72300000, 14.74900000);
  mol.emplace_back(AtomicNumber(7 ), 22.37500000,  40.72000000, 14.90700000);
  mol.emplace_back(AtomicNumber(7 ), 29.42600000,  38.43000000, 15.44600000);
  mol.emplace_back(AtomicNumber(6 ), 30.22500000,  38.64300000, 16.66200000);
  mol.emplace_back(AtomicNumber(6 ), 29.66400000,  39.83900000, 17.43400000);
  mol.emplace_back(AtomicNumber(8 ), 28.85000000,  40.56500000, 16.85900000);
  mol.emplace_back(AtomicNumber(6 ), 31.74400000,  38.87900000, 16.29900000);
  mol.emplace_back(AtomicNumber(8 ), 31.73700000,  40.25700000, 15.82400000);
  mol.emplace_back(AtomicNumber(6 ), 32.26000000,  37.96900000, 15.17100000);
  mol.emplace_back(AtomicNumber(7 ), 30.13200000,  40.06900000, 18.64200000);
  mol.emplace_back(AtomicNumber(6 ), 29.60700000,  41.18000000, 19.46700000);
  mol.emplace_back(AtomicNumber(6 ), 30.07500000,  42.53800000, 18.98400000);
  mol.emplace_back(AtomicNumber(8 ), 29.58600000,  43.57000000, 19.48300000);
  mol.emplace_back(AtomicNumber(6 ), 29.91900000,  40.89000000, 20.93800000);
  mol.emplace_back(AtomicNumber(6 ), 29.18300000,  39.72200000, 21.58100000);
  mol.emplace_back(AtomicNumber(6 ), 29.30800000,  39.75000000, 23.09500000);
  mol.emplace_back(AtomicNumber(6 ), 27.70000000,  39.72100000, 21.22800000);
  mol.emplace_back(AtomicNumber(7 ), 30.99100000,  42.57100000, 17.99800000);
  mol.emplace_back(AtomicNumber(6 ), 31.42200000,  43.94000000, 17.55300000);
  mol.emplace_back(AtomicNumber(6 ), 30.75500000,  44.35100000, 16.27700000);
  mol.emplace_back(AtomicNumber(8 ), 31.20700000,  45.26800000, 15.56600000);
  mol.emplace_back(AtomicNumber(6 ), 32.97900000,  43.91800000, 17.44500000);
  mol.emplace_back(AtomicNumber(8 ), 33.17400000,  43.06700000, 16.26500000);
  mol.emplace_back(AtomicNumber(6 ), 33.65700000,  43.31900000, 18.67200000);
  mol.emplace_back(AtomicNumber(7 ), 29.72100000,  43.67300000, 15.88500000);
  mol.emplace_back(AtomicNumber(6 ), 28.97800000,  43.96000000, 14.67800000);
  mol.emplace_back(AtomicNumber(6 ), 29.60400000,  43.50700000, 13.39300000);
  mol.emplace_back(AtomicNumber(8 ), 29.21900000,  43.98100000, 12.30100000);
  mol.emplace_back(AtomicNumber(7 ), 30.56300000,  42.62300000, 13.49500000);
  mol.emplace_back(AtomicNumber(6 ), 31.19100000,  42.01200000, 12.33100000);
  mol.emplace_back(AtomicNumber(6 ), 30.45900000,  40.66600000, 12.13000000);
  mol.emplace_back(AtomicNumber(8 ), 30.25300000,  39.99100000, 13.13300000);
  mol.emplace_back(AtomicNumber(6 ), 32.67200000,  41.71700000, 12.50500000);
  mol.emplace_back(AtomicNumber(6 ), 33.28000000,  41.08600000, 11.22700000);
  mol.emplace_back(AtomicNumber(6 ), 34.76200000,  40.79900000, 11.47000000);
  mol.emplace_back(AtomicNumber(6 ), 35.61400000,  40.84700000, 10.24000000);
  mol.emplace_back(AtomicNumber(7 ), 35.10000000,  40.07300000,  9.10100000);
  mol.emplace_back(AtomicNumber(7 ), 30.16300000,  40.33800000, 10.88600000);
  mol.emplace_back(AtomicNumber(6 ), 29.54200000,  39.02000000, 10.65300000);
  mol.emplace_back(AtomicNumber(6 ), 30.49400000,  38.26100000,  9.72900000);
  mol.emplace_back(AtomicNumber(8 ), 30.84900000,  38.85000000,  8.70600000);
  mol.emplace_back(AtomicNumber(6 ), 28.11300000,  39.04900000, 10.01500000);
  mol.emplace_back(AtomicNumber(8 ), 27.28000000,  39.72200000, 10.99600000);
  mol.emplace_back(AtomicNumber(6 ), 27.58800000,  37.63500000,  9.71500000);
  mol.emplace_back(AtomicNumber(7 ), 30.79500000,  37.01500000, 10.09500000);
  mol.emplace_back(AtomicNumber(6 ), 31.72000000,  36.28900000,  9.17600000);
  mol.emplace_back(AtomicNumber(6 ), 30.95500000,  35.21100000,  8.45900000);
  mol.emplace_back(AtomicNumber(8 ), 30.02500000,  34.61800000,  9.04000000);
  mol.emplace_back(AtomicNumber(6 ), 32.99500000,  35.88300000,  9.93400000);
  mol.emplace_back(AtomicNumber(6 ), 33.30600000,  34.38100000,  9.84000000);
  mol.emplace_back(AtomicNumber(6 ), 33.10900000,  36.38100000, 11.43500000);
  mol.emplace_back(AtomicNumber(6 ), 34.53500000,  34.02800000, 10.72000000);
  mol.emplace_back(AtomicNumber(7 ), 31.24400000,  34.98600000,  7.19700000);
  mol.emplace_back(AtomicNumber(6 ), 30.50500000,  33.88400000,  6.51200000);
  mol.emplace_back(AtomicNumber(6 ), 31.40900000,  32.68000000,  6.44600000);
  mol.emplace_back(AtomicNumber(8 ), 32.61900000,  32.81200000,  6.12500000);
  mol.emplace_back(AtomicNumber(6 ), 30.09100000,  34.39300000,  5.07800000);
  mol.emplace_back(AtomicNumber(8 ), 31.44000000,  34.51300000,  4.48700000);
  mol.emplace_back(AtomicNumber(6 ), 29.42000000,  35.75600000,  5.11900000);
  mol.emplace_back(AtomicNumber(7 ), 30.88400000,  31.48500000,  6.66600000);
  mol.emplace_back(AtomicNumber(6 ), 31.67700000,  30.27500000,  6.63900000);
  mol.emplace_back(AtomicNumber(6 ), 31.02200000,  29.28800000,  5.66500000);
  mol.emplace_back(AtomicNumber(8 ), 29.80900000,  29.39500000,  5.54500000);
  mol.emplace_back(AtomicNumber(6 ), 31.56200000,  29.68600000,  8.04500000);
  mol.emplace_back(AtomicNumber(6 ), 32.63100000,  29.44400000,  9.06000000);
  mol.emplace_back(AtomicNumber(6 ), 33.81400000,  30.39000000,  9.03000000);
  mol.emplace_back(AtomicNumber(6 ), 31.94500000,  29.44900000, 10.43600000);
  mol.emplace_back(AtomicNumber(7 ), 31.83400000,  28.41200000,  5.12500000);
  mol.emplace_back(AtomicNumber(6 ), 31.22000000,  27.34100000,  4.27500000);
  mol.emplace_back(AtomicNumber(6 ), 31.44000000,  26.07900000,  5.08000000);
  mol.emplace_back(AtomicNumber(8 ), 32.57600000,  25.80200000,  5.46100000);
  mol.emplace_back(AtomicNumber(6 ), 31.82700000,  27.26200000,  2.89400000);
  mol.emplace_back(AtomicNumber(6 ), 31.36300000,  28.41000000,  1.96200000);
  mol.emplace_back(AtomicNumber(6 ), 31.67100000,  28.29100000,  0.49800000);
  mol.emplace_back(AtomicNumber(8 ), 30.86900000,  28.62100000, -0.36600000);
  mol.emplace_back(AtomicNumber(8 ), 32.83500000,  27.86100000,  0.27800000);
  mol.emplace_back(AtomicNumber(7 ), 30.31000000,  25.45800000,  5.38400000);
  mol.emplace_back(AtomicNumber(6 ), 30.28800000,  24.24500000,  6.19300000);
  mol.emplace_back(AtomicNumber(6 ), 29.27900000,  23.22700000,  5.64100000);
  mol.emplace_back(AtomicNumber(8 ), 28.47800000,  23.52200000,  4.72500000);
  mol.emplace_back(AtomicNumber(6 ), 29.90300000,  24.59000000,  7.66500000);
  mol.emplace_back(AtomicNumber(6 ), 30.86200000,  25.49600000,  8.38900000);
  mol.emplace_back(AtomicNumber(6 ), 28.47600000,  25.13500000,  7.70500000);
  mol.emplace_back(AtomicNumber(7 ), 29.38000000,  22.05700000,  6.23200000);
  mol.emplace_back(AtomicNumber(6 ), 28.46800000,  20.94000000,  5.98000000);
  mol.emplace_back(AtomicNumber(6 ), 27.81900000,  20.60900000,  7.31600000);
  mol.emplace_back(AtomicNumber(8 ), 28.44900000,  20.67400000,  8.36000000);
  mol.emplace_back(AtomicNumber(6 ), 29.21300000,  19.69700000,  5.50600000);
  mol.emplace_back(AtomicNumber(6 ), 29.72800000,  19.75500000,  4.06000000);
  mol.emplace_back(AtomicNumber(6 ), 28.75400000,  20.06100000,  2.97800000);
  mol.emplace_back(AtomicNumber(8 ), 27.54600000,  19.99200000,  2.98500000);
  mol.emplace_back(AtomicNumber(8 ), 29.33600000,  20.42300000,  1.90400000);
  mol.emplace_back(AtomicNumber(7 ), 26.55900000,  20.22000000,  7.28800000);
  mol.emplace_back(AtomicNumber(6 ), 25.82900000,  19.82500000,  8.49400000);
  mol.emplace_back(AtomicNumber(6 ), 26.54100000,  18.73200000,  9.25100000);
  mol.emplace_back(AtomicNumber(8 ), 26.33300000,  18.53600000, 10.45700000);
  mol.emplace_back(AtomicNumber(6 ), 24.46900000,  19.33200000,  7.95200000);
  mol.emplace_back(AtomicNumber(6 ), 24.29900000,  20.13400000,  6.70400000);
  mol.emplace_back(AtomicNumber(6 ), 25.71400000,  20.10800000,  6.07300000);
  mol.emplace_back(AtomicNumber(7 ), 27.36100000,  17.95900000,  8.55900000);
  mol.emplace_back(AtomicNumber(6 ), 28.05400000,  16.83500000,  9.21000000);
  mol.emplace_back(AtomicNumber(6 ), 29.25800000,  17.31800000,  9.98400000);
  mol.emplace_back(AtomicNumber(8 ), 29.93000000,  16.47700000, 10.60600000);
  mol.emplace_back(AtomicNumber(6 ), 28.52300000,  15.82000000,  8.18200000);
  mol.emplace_back(AtomicNumber(8 ), 28.94600000,  16.44500000,  6.96700000);
  mol.emplace_back(AtomicNumber(7 ), 29.59900000,  18.59900000,  9.82800000);
  mol.emplace_back(AtomicNumber(6 ), 30.79600000,  19.08300000, 10.56600000);
  mol.emplace_back(AtomicNumber(6 ), 30.49100000,  19.16200000, 12.04000000);
  mol.emplace_back(AtomicNumber(8 ), 29.36700000,  19.52300000, 12.44100000);
  mol.emplace_back(AtomicNumber(6 ), 31.15500000,  20.51500000, 10.04800000);
  mol.emplace_back(AtomicNumber(6 ), 31.92300000,  20.43600000,  8.75500000);
  mol.emplace_back(AtomicNumber(8 ), 32.49300000,  19.37400000,  8.45600000);
  mol.emplace_back(AtomicNumber(8 ), 31.83800000,  21.40200000,  7.96800000);
  mol.emplace_back(AtomicNumber(7 ), 31.51000000,  18.93600000, 12.85200000);
  mol.emplace_back(AtomicNumber(6 ), 31.39800000,  19.06400000, 14.28600000);
  mol.emplace_back(AtomicNumber(6 ), 31.59300000,  20.55300000, 14.65500000);
  mol.emplace_back(AtomicNumber(8 ), 32.15900000,  21.31100000, 13.86100000);
  mol.emplace_back(AtomicNumber(6 ), 32.49200000,  18.19300000, 14.99500000);
  mol.emplace_back(AtomicNumber(8 ), 33.77800000,  18.73900000, 14.51600000);
  mol.emplace_back(AtomicNumber(6 ), 32.35200000,  16.70000000, 14.63000000);
  mol.emplace_back(AtomicNumber(7 ), 31.11300000,  20.86300000, 15.86000000);
  mol.emplace_back(AtomicNumber(6 ), 31.28800000,  22.20100000, 16.41700000);
  mol.emplace_back(AtomicNumber(6 ), 32.77600000,  22.51900000, 16.57700000);
  mol.emplace_back(AtomicNumber(8 ), 33.23300000,  23.65900000, 16.38400000);
  mol.emplace_back(AtomicNumber(6 ), 30.52000000,  22.30000000, 17.76400000);
  mol.emplace_back(AtomicNumber(6 ), 29.00600000,  22.04300000, 17.44200000);
  mol.emplace_back(AtomicNumber(6 ), 30.83200000,  23.69900000, 18.35800000);
  mol.emplace_back(AtomicNumber(6 ), 28.40700000,  22.94800000, 16.36600000);
  mol.emplace_back(AtomicNumber(7 ), 33.54800000,  21.52600000, 16.95000000);
  mol.emplace_back(AtomicNumber(6 ), 35.03100000,  21.72200000, 17.06900000);
  mol.emplace_back(AtomicNumber(6 ), 35.61500000,  22.19000000, 15.75900000);
  mol.emplace_back(AtomicNumber(8 ), 36.53200000,  23.04600000, 15.72400000);
  mol.emplace_back(AtomicNumber(6 ), 35.66700000,  20.38300000, 17.44700000);
  mol.emplace_back(AtomicNumber(6 ), 37.12800000,  20.29300000, 17.87200000);
  mol.emplace_back(AtomicNumber(6 ), 37.56100000,  18.85100000, 18.08200000);
  mol.emplace_back(AtomicNumber(8 ), 37.75800000,  18.02400000, 17.19500000);
  mol.emplace_back(AtomicNumber(8 ), 37.62800000,  18.59900000, 19.31300000);
  mol.emplace_back(AtomicNumber(7 ), 35.13900000,  21.62400000, 14.66200000);
  mol.emplace_back(AtomicNumber(6 ), 35.59000000,  21.94500000, 13.30200000);
  mol.emplace_back(AtomicNumber(6 ), 35.23800000,  23.38200000, 12.92000000);
  mol.emplace_back(AtomicNumber(8 ), 36.06600000,  24.10900000, 12.33300000);
  mol.emplace_back(AtomicNumber(6 ), 35.06400000,  20.95700000, 12.25500000);
  mol.emplace_back(AtomicNumber(6 ), 35.54100000,  21.41800000, 10.87100000);
  mol.emplace_back(AtomicNumber(8 ), 36.77200000,  21.62300000, 10.67600000);
  mol.emplace_back(AtomicNumber(7 ), 34.62800000,  21.59500000,  9.92000000);
  mol.emplace_back(AtomicNumber(7 ), 34.00700000,  23.74500000, 13.25000000);
  mol.emplace_back(AtomicNumber(6 ), 33.53300000,  25.09700000, 12.97800000);
  mol.emplace_back(AtomicNumber(6 ), 34.44100000,  26.09900000, 13.68400000);
  mol.emplace_back(AtomicNumber(8 ), 34.88300000,  27.09000000, 13.09300000);
  mol.emplace_back(AtomicNumber(6 ), 32.06000000,  25.25700000, 13.36400000);
  mol.emplace_back(AtomicNumber(6 ), 31.68400000,  26.74900000, 13.34200000);
  mol.emplace_back(AtomicNumber(6 ), 31.15200000,  24.42100000, 12.47700000);
  mol.emplace_back(AtomicNumber(7 ), 34.73400000,  25.82200000, 14.94900000);
  mol.emplace_back(AtomicNumber(6 ), 35.59600000,  26.71500000, 15.73600000);
  mol.emplace_back(AtomicNumber(6 ), 36.97500000,  26.82600000, 15.10700000);
  mol.emplace_back(AtomicNumber(8 ), 37.57900000,  27.92600000, 15.15900000);
  mol.emplace_back(AtomicNumber(6 ), 35.71500000,  26.20300000, 17.17200000);
  mol.emplace_back(AtomicNumber(6 ), 34.34300000,  26.44500000, 17.89800000);
  mol.emplace_back(AtomicNumber(6 ), 34.50900000,  26.07700000, 19.36000000);
  mol.emplace_back(AtomicNumber(6 ), 33.20600000,  26.31100000, 20.12200000);
  mol.emplace_back(AtomicNumber(7 ), 33.45500000,  25.91000000, 21.54600000);
  mol.emplace_back(AtomicNumber(7 ), 37.49900000,  25.74300000, 14.57100000);
  mol.emplace_back(AtomicNumber(6 ), 38.79400000,  25.76100000, 13.88000000);
  mol.emplace_back(AtomicNumber(6 ), 38.72800000,  26.59100000, 12.61100000);
  mol.emplace_back(AtomicNumber(8 ), 39.70400000,  27.34600000, 12.27700000);
  mol.emplace_back(AtomicNumber(6 ), 39.28500000,  24.33600000, 13.56600000);
  mol.emplace_back(AtomicNumber(7 ), 37.63300000,  26.54300000, 11.86700000);
  mol.emplace_back(AtomicNumber(6 ), 37.47100000,  27.39100000, 10.66800000);
  mol.emplace_back(AtomicNumber(6 ), 37.44100000,  28.88200000, 11.05200000);
  mol.emplace_back(AtomicNumber(8 ), 38.02000000,  29.77200000, 10.38200000);
  mol.emplace_back(AtomicNumber(6 ), 36.19300000,  27.05800000,  9.91100000);
  mol.emplace_back(AtomicNumber(6 ), 36.15300000,  25.62000000,  9.40900000);
  mol.emplace_back(AtomicNumber(6 ), 34.75800000,  25.28000000,  8.90000000);
  mol.emplace_back(AtomicNumber(6 ), 34.79300000,  24.26400000,  7.76700000);
  mol.emplace_back(AtomicNumber(7 ), 34.91400000,  24.94400000,  6.44100000);
  mol.emplace_back(AtomicNumber(7 ), 36.81100000,  29.17000000, 12.19200000);
  mol.emplace_back(AtomicNumber(6 ), 36.73100000,  30.57000000, 12.64500000);
  mol.emplace_back(AtomicNumber(6 ), 38.14800000,  30.98100000, 13.06900000);
  mol.emplace_back(AtomicNumber(8 ), 38.54400000,  32.15000000, 12.85600000);
  mol.emplace_back(AtomicNumber(6 ), 35.70800000,  30.77600000, 13.80600000);
  mol.emplace_back(AtomicNumber(6 ), 34.22800000,  30.63000000, 13.31900000);
  mol.emplace_back(AtomicNumber(6 ), 35.87400000,  32.13800000, 14.51200000);
  mol.emplace_back(AtomicNumber(6 ), 33.28400000,  30.50400000, 14.55200000);
  mol.emplace_back(AtomicNumber(7 ), 38.88300000,  30.11000000, 13.71300000);
  mol.emplace_back(AtomicNumber(6 ), 40.26900000,  30.50800000, 14.11500000);
  mol.emplace_back(AtomicNumber(6 ), 41.09200000,  30.80800000, 12.85100000);
  mol.emplace_back(AtomicNumber(8 ), 41.82800000,  31.80800000, 12.68100000);
  mol.emplace_back(AtomicNumber(6 ), 40.99600000,  29.39900000, 14.86500000);
  mol.emplace_back(AtomicNumber(6 ), 42.44500000,  29.84800000, 15.18200000);
  mol.emplace_back(AtomicNumber(6 ), 43.09000000,  28.82800000, 16.09500000);
  mol.emplace_back(AtomicNumber(8 ), 42.77000000,  27.65500000, 15.90600000);
  mol.emplace_back(AtomicNumber(7 ), 43.89800000,  29.25200000, 17.05000000);
  mol.emplace_back(AtomicNumber(7 ), 41.00100000,  29.87800000, 11.93100000);
  mol.emplace_back(AtomicNumber(6 ), 41.71800000,  30.02200000, 10.64300000);
  mol.emplace_back(AtomicNumber(6 ), 41.39900000,  31.33800000,  9.96700000);
  mol.emplace_back(AtomicNumber(8 ), 42.26000000,  32.03600000,  9.38100000);
  mol.emplace_back(AtomicNumber(6 ), 41.39800000,  28.78000000,  9.81000000);
  mol.emplace_back(AtomicNumber(6 ), 42.62600000,  28.55700000,  8.92800000);
  mol.emplace_back(AtomicNumber(8 ), 43.66600000,  28.26200000,  9.53900000);
  mol.emplace_back(AtomicNumber(8 ), 42.43000000,  28.81200000,  7.72800000);
  mol.emplace_back(AtomicNumber(7 ), 40.11700000,  31.75000000,  9.98800000);
  mol.emplace_back(AtomicNumber(6 ), 39.80800000,  32.99400000,  9.23300000);
  mol.emplace_back(AtomicNumber(6 ), 39.83700000,  34.27100000,  9.99500000);
  mol.emplace_back(AtomicNumber(8 ), 40.16400000,  35.32300000,  9.34500000);
  mol.emplace_back(AtomicNumber(6 ), 38.61500000,  32.80100000,  8.32000000);
  mol.emplace_back(AtomicNumber(6 ), 37.22000000,  32.82200000,  8.82700000);
  mol.emplace_back(AtomicNumber(6 ), 36.35100000,  33.61300000,  7.83800000);
  mol.emplace_back(AtomicNumber(6 ), 36.32200000,  32.94400000,  6.47700000);
  mol.emplace_back(AtomicNumber(7 ), 35.76800000,  33.94500000,  5.48900000);
  mol.emplace_back(AtomicNumber(7 ), 39.65500000,  34.33500000, 11.28500000);
  mol.emplace_back(AtomicNumber(6 ), 39.67600000,  35.54700000, 12.07200000);
  mol.emplace_back(AtomicNumber(6 ), 40.67500000,  35.52700000, 13.20000000);
  mol.emplace_back(AtomicNumber(8 ), 40.81400000,  36.52800000, 13.91100000);
  mol.emplace_back(AtomicNumber(6 ), 38.29000000,  35.81400000, 12.69800000);
  mol.emplace_back(AtomicNumber(6 ), 37.15600000,  35.98500000, 11.68800000);
  mol.emplace_back(AtomicNumber(6 ), 37.19200000,  37.36100000, 11.03300000);
  mol.emplace_back(AtomicNumber(8 ), 37.51900000,  38.36000000, 11.64500000);
  mol.emplace_back(AtomicNumber(8 ), 36.86100000,  37.32000000,  9.82200000);
  mol.emplace_back(AtomicNumber(7 ), 41.31700000,  34.39300000, 13.43200000);
  mol.emplace_back(AtomicNumber(6 ), 42.34500000,  34.26900000, 14.43100000);
  mol.emplace_back(AtomicNumber(6 ), 41.94900000,  34.07600000, 15.84200000);
  mol.emplace_back(AtomicNumber(8 ), 42.82900000,  34.00000000, 16.73900000);
  mol.emplace_back(AtomicNumber(7 ), 40.64200000,  33.91600000, 16.11200000);
  mol.emplace_back(AtomicNumber(6 ), 40.22600000,  33.71600000, 17.50900000);
  mol.emplace_back(AtomicNumber(6 ), 40.44900000,  32.27800000, 17.94500000);
  mol.emplace_back(AtomicNumber(8 ), 39.93600000,  31.33600000, 17.31500000);
  mol.emplace_back(AtomicNumber(6 ), 38.69300000,  34.10600000, 17.59500000);
  mol.emplace_back(AtomicNumber(6 ), 38.47100000,  35.54600000, 17.04500000);
  mol.emplace_back(AtomicNumber(6 ), 38.14600000,  33.93200000, 19.02700000);
  mol.emplace_back(AtomicNumber(6 ), 36.95800000,  35.74600000, 16.68000000);
  mol.emplace_back(AtomicNumber(7 ), 41.18900000,  32.08500000, 19.03100000);
  mol.emplace_back(AtomicNumber(6 ), 41.46100000,  30.75100000, 19.59400000);
  mol.emplace_back(AtomicNumber(6 ), 40.16800000,  30.02600000, 19.91800000);
  mol.emplace_back(AtomicNumber(8 ), 39.26400000,  30.66200000, 20.52100000);
  mol.emplace_back(AtomicNumber(6 ), 42.19500000,  31.14200000, 20.91300000);
  mol.emplace_back(AtomicNumber(6 ), 42.90400000,  32.41400000, 20.55300000);
  mol.emplace_back(AtomicNumber(6 ), 41.82200000,  33.18800000, 19.81300000);
  mol.emplace_back(AtomicNumber(7 ), 40.05900000,  28.75800000, 19.60700000);
  mol.emplace_back(AtomicNumber(6 ), 38.81700000,  28.02000000, 19.88900000);
  mol.emplace_back(AtomicNumber(6 ), 38.42100000,  28.04800000, 21.34100000);
  mol.emplace_back(AtomicNumber(8 ), 37.21300000,  28.03600000, 21.70400000);
  mol.emplace_back(AtomicNumber(6 ), 39.09000000,  26.62900000, 19.32500000);
  mol.emplace_back(AtomicNumber(6 ), 40.08200000,  26.90400000, 18.19800000);
  mol.emplace_back(AtomicNumber(6 ), 41.03500000,  27.90900000, 18.87900000);
  mol.emplace_back(AtomicNumber(7 ), 39.37400000,  28.09000000, 22.24000000);
  mol.emplace_back(AtomicNumber(6 ), 39.06300000,  28.06300000, 23.69500000);
  mol.emplace_back(AtomicNumber(6 ), 38.36500000,  29.33500000, 24.15900000);
  mol.emplace_back(AtomicNumber(8 ), 37.68400000,  29.39000000, 25.22100000);
  mol.emplace_back(AtomicNumber(6 ), 40.34000000,  27.69200000, 24.46800000);
  mol.emplace_back(AtomicNumber(6 ), 40.55900000,  28.58500000, 25.67500000);
  mol.emplace_back(AtomicNumber(8 ), 40.71600000,  29.80900000, 25.45600000);
  mol.emplace_back(AtomicNumber(8 ), 40.54900000,  28.09000000, 26.84000000);
  mol.emplace_back(AtomicNumber(7 ), 38.41900000,  30.37300000, 23.34100000);
  mol.emplace_back(AtomicNumber(6 ), 37.73800000,  31.63700000, 23.71200000);
  mol.emplace_back(AtomicNumber(6 ), 36.33400000,  31.74200000, 23.08700000);
  mol.emplace_back(AtomicNumber(8 ), 35.57400000,  32.61800000, 23.48300000);
  mol.emplace_back(AtomicNumber(6 ), 38.52800000,  32.85400000, 23.18200000);
  mol.emplace_back(AtomicNumber(6 ), 39.91900000,  32.85400000, 23.84000000);
  mol.emplace_back(AtomicNumber(6 ), 40.76000000,  34.03600000, 23.39400000);
  mol.emplace_back(AtomicNumber(8 ), 41.97500000,  34.00800000, 23.62400000);
  mol.emplace_back(AtomicNumber(7 ), 40.14000000,  35.00700000, 22.77500000);
  mol.emplace_back(AtomicNumber(7 ), 36.00000000,  30.86000000, 22.17200000);
  mol.emplace_back(AtomicNumber(6 ), 34.73800000,  30.87500000, 21.47300000);
  mol.emplace_back(AtomicNumber(6 ), 33.58900000,  30.18900000, 22.18100000);
  mol.emplace_back(AtomicNumber(8 ), 33.58000000,  29.00900000, 22.49900000);
  mol.emplace_back(AtomicNumber(6 ), 34.87600000,  30.23700000, 20.06600000);
  mol.emplace_back(AtomicNumber(6 ), 36.01200000,  30.86000000, 19.22100000);
  mol.emplace_back(AtomicNumber(6 ), 36.08300000,  30.19400000, 17.87500000);
  mol.emplace_back(AtomicNumber(8 ), 35.04800000,  29.70200000, 17.39300000);
  mol.emplace_back(AtomicNumber(7 ), 37.22800000,  30.12600000, 17.23300000);
  mol.emplace_back(AtomicNumber(7 ), 32.47800000,  30.91700000, 22.26900000);
  mol.emplace_back(AtomicNumber(6 ), 31.20000000,  30.32900000, 22.78000000);
  mol.emplace_back(AtomicNumber(6 ), 30.21000000,  30.50900000, 21.65000000);
  mol.emplace_back(AtomicNumber(8 ), 29.97800000,  31.72600000, 21.26900000);
  mol.emplace_back(AtomicNumber(6 ), 30.84700000,  30.93100000, 24.11800000);
  mol.emplace_back(AtomicNumber(6 ), 29.41200000,  30.79600000, 24.59800000);
  mol.emplace_back(AtomicNumber(6 ), 29.27100000,  31.31400000, 26.01600000);
  mol.emplace_back(AtomicNumber(7 ), 27.87500000,  31.31700000, 26.44300000);
  mol.emplace_back(AtomicNumber(6 ), 27.13200000,  32.42300000, 26.57400000);
  mol.emplace_back(AtomicNumber(7 ), 27.63000000,  33.65600000, 26.46100000);
  mol.emplace_back(AtomicNumber(7 ), 25.81000000,  32.29900000, 26.73200000);
  mol.emplace_back(AtomicNumber(7 ), 29.69400000,  29.43600000, 21.05400000);
  mol.emplace_back(AtomicNumber(6 ), 28.76200000,  29.57300000, 19.90600000);
  mol.emplace_back(AtomicNumber(6 ), 27.33100000,  29.31700000, 20.36400000);
  mol.emplace_back(AtomicNumber(8 ), 27.10100000,  28.34600000, 21.09700000);
  mol.emplace_back(AtomicNumber(6 ), 29.15100000,  28.65500000, 18.75500000);
  mol.emplace_back(AtomicNumber(6 ), 30.41600000,  28.91200000, 17.98000000);
  mol.emplace_back(AtomicNumber(6 ), 30.73800000,  27.69300000, 17.12200000);
  mol.emplace_back(AtomicNumber(6 ), 30.20500000,  30.16800000, 17.12900000);
  mol.emplace_back(AtomicNumber(7 ), 26.43600000,  30.23200000, 20.00400000);
  mol.emplace_back(AtomicNumber(6 ), 25.03400000,  30.17000000, 20.40100000);
  mol.emplace_back(AtomicNumber(6 ), 24.10100000,  30.14900000, 19.19600000);
  mol.emplace_back(AtomicNumber(8 ), 24.19600000,  30.94800000, 18.28700000);
  mol.emplace_back(AtomicNumber(6 ), 24.63900000,  31.42600000, 21.28600000);
  mol.emplace_back(AtomicNumber(6 ), 25.64600000,  31.67000000, 22.42100000);
  mol.emplace_back(AtomicNumber(6 ), 23.18100000,  31.30900000, 21.82400000);
  mol.emplace_back(AtomicNumber(6 ), 25.77800000,  30.43600000, 23.35600000);
  mol.emplace_back(AtomicNumber(7 ), 23.14100000,  29.18700000, 19.24100000);
  mol.emplace_back(AtomicNumber(6 ), 22.12600000,  29.06200000, 18.18300000);
  mol.emplace_back(AtomicNumber(6 ), 20.83500000,  28.62900000, 18.90400000);
  mol.emplace_back(AtomicNumber(8 ), 20.82100000,  27.73400000, 19.74900000);
  mol.emplace_back(AtomicNumber(6 ), 22.49400000,  28.05700000, 17.10900000);
  mol.emplace_back(AtomicNumber(6 ), 21.44700000,  27.86900000, 16.02600000);
  mol.emplace_back(AtomicNumber(6 ), 21.32500000,  28.81300000, 15.00500000);
  mol.emplace_back(AtomicNumber(6 ), 20.63800000,  26.73500000, 16.05300000);
  mol.emplace_back(AtomicNumber(6 ), 20.36900000,  28.64800000, 14.00100000);
  mol.emplace_back(AtomicNumber(6 ), 19.67700000,  26.53900000, 15.05100000);
  mol.emplace_back(AtomicNumber(6 ), 19.59300000,  27.46500000, 14.02100000);
  mol.emplace_back(AtomicNumber(7 ), 19.81000000,  29.37800000, 18.57800000);
  mol.emplace_back(AtomicNumber(6 ), 18.44300000,  29.14300000, 19.08300000);
  mol.emplace_back(AtomicNumber(6 ), 18.45300000,  28.94100000, 20.59100000);
  mol.emplace_back(AtomicNumber(8 ), 17.86000000,  27.99400000, 21.12800000);
  mol.emplace_back(AtomicNumber(6 ), 17.86400000,  27.97700000, 18.34600000);
  mol.emplace_back(AtomicNumber(7 ), 19.17200000,  29.80800000, 21.24300000);
  mol.emplace_back(AtomicNumber(6 ), 19.39900000,  29.89400000, 22.65500000);
  mol.emplace_back(AtomicNumber(6 ), 20.08300000,  28.72900000, 23.32100000);
  mol.emplace_back(AtomicNumber(8 ), 19.99100000,  28.58400000, 24.56100000);
  mol.emplace_back(AtomicNumber(7 ), 20.80100000,  27.93100000, 22.57800000);
  mol.emplace_back(AtomicNumber(6 ), 21.55000000,  26.79600000, 23.13300000);
  mol.emplace_back(AtomicNumber(6 ), 23.04600000,  27.08700000, 22.91300000);
  mol.emplace_back(AtomicNumber(8 ), 23.38300000,  27.62700000, 21.87000000);
  mol.emplace_back(AtomicNumber(6 ), 21.24200000,  25.51900000, 22.39100000);
  mol.emplace_back(AtomicNumber(6 ), 19.76200000,  25.07700000, 22.45500000);
  mol.emplace_back(AtomicNumber(6 ), 19.63400000,  23.88500000, 21.53100000);
  mol.emplace_back(AtomicNumber(6 ), 18.79100000,  24.22100000, 20.31300000);
  mol.emplace_back(AtomicNumber(7 ), 17.44000000,  24.65500000, 20.82700000);
  mol.emplace_back(AtomicNumber(7 ), 23.88000000,  26.72700000, 23.85100000);
  mol.emplace_back(AtomicNumber(6 ), 25.34900000,  26.87200000, 23.64300000);
  mol.emplace_back(AtomicNumber(6 ), 25.74300000,  25.58600000, 22.92200000);
  mol.emplace_back(AtomicNumber(8 ), 25.32500000,  24.48900000, 23.37800000);
  mol.emplace_back(AtomicNumber(6 ), 26.07000000,  27.02500000, 24.96000000);
  mol.emplace_back(AtomicNumber(6 ), 27.55300000,  27.35600000, 24.69500000);
  mol.emplace_back(AtomicNumber(6 ), 28.26200000,  27.57600000, 26.02000000);
  mol.emplace_back(AtomicNumber(8 ), 29.18900000,  26.84000000, 26.33500000);
  mol.emplace_back(AtomicNumber(7 ), 27.77700000,  28.58500000, 26.73900000);
  mol.emplace_back(AtomicNumber(7 ), 26.46500000,  25.68900000, 21.83300000);
  mol.emplace_back(AtomicNumber(6 ), 26.82600000,  24.52100000, 21.01200000);
  mol.emplace_back(AtomicNumber(6 ), 27.99400000,  23.78100000, 21.64300000);
  mol.emplace_back(AtomicNumber(8 ), 28.90400000,  24.44400000, 22.09800000);
  mol.emplace_back(AtomicNumber(6 ), 27.04300000,  24.99200000, 19.57100000);
  mol.emplace_back(AtomicNumber(6 ), 25.93100000,  25.84400000, 18.95900000);
  mol.emplace_back(AtomicNumber(6 ), 26.20300000,  26.08300000, 17.47100000);
  mol.emplace_back(AtomicNumber(6 ), 24.57700000,  25.19000000, 19.07900000);
  mol.emplace_back(AtomicNumber(7 ), 27.94200000,  22.44800000, 21.64800000);
  mol.emplace_back(AtomicNumber(6 ), 29.01500000,  21.65700000, 22.28800000);
  mol.emplace_back(AtomicNumber(6 ), 29.94200000,  21.10600000, 21.24000000);
  mol.emplace_back(AtomicNumber(8 ), 29.47000000,  20.67700000, 20.19000000);
  mol.emplace_back(AtomicNumber(6 ), 28.34800000,  20.54000000, 23.06600000);
  mol.emplace_back(AtomicNumber(6 ), 29.24700000,  19.45600000, 23.70500000);
  mol.emplace_back(AtomicNumber(6 ), 28.72200000,  19.04700000, 25.06600000);
  mol.emplace_back(AtomicNumber(8 ), 29.13900000,  18.13200000, 25.74600000);
  mol.emplace_back(AtomicNumber(8 ), 27.77700000,  19.84200000, 25.36700000);
  mol.emplace_back(AtomicNumber(7 ), 31.23300000,  21.09000000, 21.45900000);
  mol.emplace_back(AtomicNumber(6 ), 32.26200000,  20.67000000, 20.51400000);
  mol.emplace_back(AtomicNumber(6 ), 32.12800000,  19.36400000, 19.75000000);
  mol.emplace_back(AtomicNumber(8 ), 32.54600000,  19.31700000, 18.55800000);
  mol.emplace_back(AtomicNumber(6 ), 33.63800000,  20.71600000, 21.24200000);
  mol.emplace_back(AtomicNumber(6 ), 34.17400000,  22.12900000, 21.35400000);
  mol.emplace_back(AtomicNumber(8 ), 35.25200000,  22.32200000, 21.95800000);
  mol.emplace_back(AtomicNumber(8 ), 33.54400000,  23.08600000, 20.88300000);
  mol.emplace_back(AtomicNumber(7 ), 31.69700000,  18.31100000, 20.40600000);
  mol.emplace_back(AtomicNumber(6 ), 31.56800000,  16.96200000, 19.82500000);
  mol.emplace_back(AtomicNumber(6 ), 30.32000000,  16.69800000, 19.05100000);
  mol.emplace_back(AtomicNumber(8 ), 30.19800000,  15.65700000, 18.36600000);
  mol.emplace_back(AtomicNumber(7 ), 29.34000000,  17.59400000, 19.07600000);
  mol.emplace_back(AtomicNumber(6 ), 28.10800000,  17.43900000, 18.27600000);
  mol.emplace_back(AtomicNumber(6 ), 28.37500000,  17.99900000, 16.88700000);
  mol.emplace_back(AtomicNumber(8 ), 29.32600000,  18.78600000, 16.69000000);
  mol.emplace_back(AtomicNumber(6 ), 26.92600000,  18.19100000, 18.89200000);
  mol.emplace_back(AtomicNumber(6 ), 26.62100000,  17.79900000, 20.35200000);
  mol.emplace_back(AtomicNumber(6 ), 26.01000000,  16.37000000, 20.28000000);
  mol.emplace_back(AtomicNumber(7 ), 26.97500000,  15.52100000, 20.94200000);
  mol.emplace_back(AtomicNumber(6 ), 27.60300000,  14.42300000, 20.65500000);
  mol.emplace_back(AtomicNumber(7 ), 27.47900000,  13.73300000, 19.53700000);
  mol.emplace_back(AtomicNumber(7 ), 28.51900000,  13.96700000, 21.55000000);
  mol.emplace_back(AtomicNumber(7 ), 27.51000000,  17.68900000, 15.95400000);
  mol.emplace_back(AtomicNumber(6 ), 27.57400000,  18.19200000, 14.56300000);
  mol.emplace_back(AtomicNumber(6 ), 26.48200000,  19.28000000, 14.43200000);
  mol.emplace_back(AtomicNumber(8 ), 25.60900000,  19.38800000, 15.28700000);
  mol.emplace_back(AtomicNumber(6 ), 27.29900000,  17.05500000, 13.53300000);
  mol.emplace_back(AtomicNumber(8 ), 25.92500000,  16.61100000, 13.91300000);
  mol.emplace_back(AtomicNumber(6 ), 28.23600000,  15.86400000, 13.55800000);
  mol.emplace_back(AtomicNumber(7 ), 26.58500000,  20.06300000, 13.37800000);
  mol.emplace_back(AtomicNumber(6 ), 25.59400000,  21.10900000, 13.07200000);
  mol.emplace_back(AtomicNumber(6 ), 24.24100000,  20.43600000, 12.85700000);
  mol.emplace_back(AtomicNumber(8 ), 23.26400000,  20.95100000, 13.32900000);
  mol.emplace_back(AtomicNumber(6 ), 26.08400000,  21.88800000, 11.83300000);
  mol.emplace_back(AtomicNumber(6 ), 27.42600000,  22.61600000, 11.90200000);
  mol.emplace_back(AtomicNumber(6 ), 27.71800000,  23.34100000, 10.57800000);
  mol.emplace_back(AtomicNumber(6 ), 27.38000000,  23.72100000, 12.95500000);
  mol.emplace_back(AtomicNumber(7 ), 24.24000000,  19.23300000, 12.24600000);
  mol.emplace_back(AtomicNumber(6 ), 22.92400000,  18.58300000, 12.02500000);
  mol.emplace_back(AtomicNumber(6 ), 22.22900000,  18.24400000, 13.32500000);
  mol.emplace_back(AtomicNumber(8 ), 20.96300000,  18.25300000, 13.39500000);
  mol.emplace_back(AtomicNumber(6 ), 23.05900000,  17.32600000, 11.15400000);
  mol.emplace_back(AtomicNumber(8 ), 23.91400000,  16.39500000, 11.75500000);
  mol.emplace_back(AtomicNumber(7 ), 22.99700000,  17.97800000, 14.36600000);
  mol.emplace_back(AtomicNumber(6 ), 22.41800000,  17.63800000, 15.69300000);
  mol.emplace_back(AtomicNumber(6 ), 21.46000000,  18.73700000, 16.16300000);
  mol.emplace_back(AtomicNumber(8 ), 20.49700000,  18.50600000, 16.90000000);
  mol.emplace_back(AtomicNumber(6 ), 23.46100000,  17.33100000, 16.74100000);
  mol.emplace_back(AtomicNumber(6 ), 24.18400000,  16.01600000, 16.61900000);
  mol.emplace_back(AtomicNumber(8 ), 25.30300000,  15.89400000, 17.15200000);
  mol.emplace_back(AtomicNumber(8 ), 23.57200000,  15.10700000, 15.97500000);
  mol.emplace_back(AtomicNumber(7 ), 21.84600000,  19.95400000, 15.90500000);
  mol.emplace_back(AtomicNumber(6 ), 21.07900000,  21.14900000, 16.25100000);
  mol.emplace_back(AtomicNumber(6 ), 20.14200000,  21.59000000, 15.14900000);
  mol.emplace_back(AtomicNumber(8 ), 19.49900000,  22.64500000, 15.32100000);
  mol.emplace_back(AtomicNumber(6 ), 22.08500000,  22.25400000, 16.58100000);
  mol.emplace_back(AtomicNumber(6 ), 22.94500000,  21.95100000, 17.78500000);
  mol.emplace_back(AtomicNumber(6 ), 24.27200000,  21.54400000, 17.64400000);
  mol.emplace_back(AtomicNumber(6 ), 22.43700000,  22.15700000, 19.06500000);
  mol.emplace_back(AtomicNumber(6 ), 25.05200000,  21.28500000, 18.77600000);
  mol.emplace_back(AtomicNumber(6 ), 23.20400000,  21.90700000, 20.19200000);
  mol.emplace_back(AtomicNumber(6 ), 24.51700000,  21.47000000, 20.03000000);
  mol.emplace_back(AtomicNumber(8 ), 25.24800000,  21.30200000, 21.19100000);
  mol.emplace_back(AtomicNumber(7 ), 19.99300000,  20.88400000, 14.04900000);
  mol.emplace_back(AtomicNumber(6 ), 19.06500000,  21.35200000, 12.99900000);
  mol.emplace_back(AtomicNumber(6 ), 19.44200000,  22.74500000, 12.51000000);
  mol.emplace_back(AtomicNumber(8 ), 18.57100000,  23.61000000, 12.28900000);
  mol.emplace_back(AtomicNumber(6 ), 17.58600000,  21.28200000, 13.46100000);
  mol.emplace_back(AtomicNumber(6 ), 16.57600000,  21.25800000, 12.31500000);
  mol.emplace_back(AtomicNumber(8 ), 15.44000000,  21.81900000, 12.37800000);
  mol.emplace_back(AtomicNumber(7 ), 16.92400000,  20.58600000, 11.21600000);
  mol.emplace_back(AtomicNumber(7 ), 20.71700000,  22.96400000, 12.26000000);
  mol.emplace_back(AtomicNumber(6 ), 21.18400000,  24.26300000, 11.69000000);
  mol.emplace_back(AtomicNumber(6 ), 21.11000000,  24.11100000, 10.17300000);
  mol.emplace_back(AtomicNumber(8 ), 21.84100000,  23.19800000,  9.68600000);
  mol.emplace_back(AtomicNumber(6 ), 22.65000000,  24.51600000, 12.17200000);
  mol.emplace_back(AtomicNumber(6 ), 22.66200000,  24.81900000, 13.69900000);
  mol.emplace_back(AtomicNumber(6 ), 23.37600000,  25.64500000, 11.40900000);
  mol.emplace_back(AtomicNumber(6 ), 24.12300000,  24.98100000, 14.19500000);
  mol.emplace_back(AtomicNumber(7 ), 20.29100000,  24.87500000,  9.50700000);
  mol.emplace_back(AtomicNumber(6 ), 20.08100000,  24.77300000,  8.03300000);
  mol.emplace_back(AtomicNumber(6 ), 20.82200000,  25.91400000,  7.33200000);
  mol.emplace_back(AtomicNumber(8 ), 21.32300000,  26.83000000,  8.00800000);
  mol.emplace_back(AtomicNumber(6 ), 18.59900000,  24.73600000,  7.72700000);
  mol.emplace_back(AtomicNumber(6 ), 17.81900000,  23.43400000,  7.90000000);
  mol.emplace_back(AtomicNumber(6 ), 16.50900000,  23.52900000,  7.11600000);
  mol.emplace_back(AtomicNumber(8 ), 15.44600000,  22.98000000,  7.43300000);
  mol.emplace_back(AtomicNumber(7 ), 16.53900000,  24.29300000,  6.00900000);
  mol.emplace_back(AtomicNumber(7 ), 20.92400000,  25.86200000,  6.00600000);
  mol.emplace_back(AtomicNumber(6 ), 21.65600000,  26.84700000,  5.24000000);
  mol.emplace_back(AtomicNumber(6 ), 21.12700000,  28.24000000,  5.57400000);
  mol.emplace_back(AtomicNumber(8 ), 19.95800000,  28.46500000,  5.84200000);
  mol.emplace_back(AtomicNumber(6 ), 21.63100000,  26.64200000,  3.73100000);
  mol.emplace_back(AtomicNumber(6 ), 20.21000000,  26.42300000,  3.17500000);
  mol.emplace_back(AtomicNumber(6 ), 20.26800000,  26.58900000,  1.65600000);
  mol.emplace_back(AtomicNumber(6 ), 19.20200000,  25.85700000,  0.89100000);
  mol.emplace_back(AtomicNumber(7 ), 17.88400000,  26.54400000,  1.07500000);
  mol.emplace_back(AtomicNumber(7 ), 22.09900000,  29.16300000,  5.60500000);
  mol.emplace_back(AtomicNumber(6 ), 21.90700000,  30.56300000,  5.88100000);
  mol.emplace_back(AtomicNumber(6 ), 21.46600000,  30.95300000,  7.26100000);
  mol.emplace_back(AtomicNumber(8 ), 21.06600000,  32.11200000,  7.53300000);
  mol.emplace_back(AtomicNumber(6 ), 21.02300000,  31.22300000,  4.78400000);
  mol.emplace_back(AtomicNumber(6 ), 21.86100000,  31.34200000,  3.47400000);
  mol.emplace_back(AtomicNumber(6 ), 21.15600000,  30.72600000,  2.31100000);
  mol.emplace_back(AtomicNumber(8 ), 19.94200000,  30.79300000,  2.17000000);
  mol.emplace_back(AtomicNumber(8 ), 21.95400000,  30.15200000,  1.53500000);
  mol.emplace_back(AtomicNumber(7 ), 21.67400000,  30.03400000,  8.19100000);
  mol.emplace_back(AtomicNumber(6 ), 21.41900000,  30.25300000,  9.62000000);
  mol.emplace_back(AtomicNumber(6 ), 22.50400000,  31.22800000, 10.13600000);
  mol.emplace_back(AtomicNumber(8 ), 23.57900000,  31.32100000,  9.55400000);
  mol.emplace_back(AtomicNumber(6 ), 21.63700000,  28.92300000, 10.35300000);
  mol.emplace_back(AtomicNumber(8 ), 20.54400000,  28.04700000, 10.05900000);
  mol.emplace_back(AtomicNumber(7 ), 22.24100000,  31.87300000, 11.24100000);
  mol.emplace_back(AtomicNumber(6 ), 23.21200000,  32.76200000, 11.89100000);
  mol.emplace_back(AtomicNumber(6 ), 23.50900000,  32.22400000, 13.29000000);
  mol.emplace_back(AtomicNumber(8 ), 22.54400000,  31.94200000, 14.03400000);
  mol.emplace_back(AtomicNumber(6 ), 22.69900000,  34.26700000, 11.98500000);
  mol.emplace_back(AtomicNumber(8 ), 22.49500000,  34.69000000, 10.58900000);
  mol.emplace_back(AtomicNumber(6 ), 23.72700000,  35.13100000, 12.72200000);
  mol.emplace_back(AtomicNumber(7 ), 24.79000000,  32.02100000, 13.61800000);
  mol.emplace_back(AtomicNumber(6 ), 25.14900000,  31.60900000, 14.98000000);
  mol.emplace_back(AtomicNumber(6 ), 25.69800000,  32.87600000, 15.66900000);
  mol.emplace_back(AtomicNumber(8 ), 26.15800000,  33.73000000, 14.89400000);
  mol.emplace_back(AtomicNumber(6 ), 26.31000000,  30.59400000, 14.96700000);
  mol.emplace_back(AtomicNumber(6 ), 26.29000000,  29.48000000, 13.96000000);
  mol.emplace_back(AtomicNumber(6 ), 27.39300000,  28.44200000, 14.22900000);
  mol.emplace_back(AtomicNumber(6 ), 24.94200000,  28.80700000, 13.95200000);
  mol.emplace_back(AtomicNumber(7 ), 25.62100000,  32.94500000, 16.95000000);
  mol.emplace_back(AtomicNumber(6 ), 26.17900000,  34.12700000, 17.65000000);
  mol.emplace_back(AtomicNumber(6 ), 27.47500000,  33.65100000, 18.30400000);
  mol.emplace_back(AtomicNumber(8 ), 27.50700000,  32.58700000, 18.95800000);
  mol.emplace_back(AtomicNumber(6 ), 25.21400000,  34.56500000, 18.78000000);
  mol.emplace_back(AtomicNumber(6 ), 23.97800000,  35.12100000, 18.12600000);
  mol.emplace_back(AtomicNumber(7 ), 23.85300000,  36.43200000, 17.78100000);
  mol.emplace_back(AtomicNumber(6 ), 22.82400000,  34.51400000, 17.78200000);
  mol.emplace_back(AtomicNumber(6 ), 22.67400000,  36.62700000, 17.20000000);
  mol.emplace_back(AtomicNumber(7 ), 22.04500000,  35.45500000, 17.17300000);
  mol.emplace_back(AtomicNumber(7 ), 28.52500000,  34.44700000, 18.18900000);
  mol.emplace_back(AtomicNumber(6 ), 29.80100000,  34.14500000, 18.82900000);
  mol.emplace_back(AtomicNumber(6 ), 30.05200000,  35.04200000, 20.00400000);
  mol.emplace_back(AtomicNumber(8 ), 30.10500000,  36.30500000, 19.78800000);
  mol.emplace_back(AtomicNumber(6 ), 30.92500000,  34.30400000, 17.75300000);
  mol.emplace_back(AtomicNumber(6 ), 32.34500000,  34.18300000, 18.35800000);
  mol.emplace_back(AtomicNumber(6 ), 32.55500000,  32.78300000, 18.87000000);
  mol.emplace_back(AtomicNumber(6 ), 33.36100000,  34.49100000, 17.24500000);
  mol.emplace_back(AtomicNumber(7 ), 30.12400000,  34.53300000, 21.19100000);
  mol.emplace_back(AtomicNumber(6 ), 30.47900000,  35.36900000, 22.37400000);
  mol.emplace_back(AtomicNumber(6 ), 31.90100000,  34.91000000, 22.72800000);
  mol.emplace_back(AtomicNumber(8 ), 32.19000000,  33.69600000, 22.63500000);
  mol.emplace_back(AtomicNumber(6 ), 29.47200000,  35.18100000, 23.49800000);
  mol.emplace_back(AtomicNumber(6 ), 29.82100000,  35.95700000, 24.76500000);
  mol.emplace_back(AtomicNumber(6 ), 28.04900000,  35.45400000, 23.07100000);
  mol.emplace_back(AtomicNumber(7 ), 32.76300000,  35.83100000, 23.09000000);
  mol.emplace_back(AtomicNumber(6 ), 34.14500000,  35.47200000, 23.48100000);
  mol.emplace_back(AtomicNumber(6 ), 34.23900000,  35.35300000, 24.97900000);
  mol.emplace_back(AtomicNumber(8 ), 33.70700000,  36.19700000, 25.72800000);
  mol.emplace_back(AtomicNumber(6 ), 35.11400000,  36.56400000, 22.90700000);
  mol.emplace_back(AtomicNumber(6 ), 35.92600000,  35.97900000, 21.73700000);
  mol.emplace_back(AtomicNumber(6 ), 35.00300000,  35.08400000, 20.92000000);
  mol.emplace_back(AtomicNumber(6 ), 36.53300000,  37.08700000, 20.91700000);
  mol.emplace_back(AtomicNumber(7 ), 34.93000000,  34.38400000, 25.45100000);
  mol.emplace_back(AtomicNumber(6 ), 35.16100000,  34.17400000, 26.89600000);
  mol.emplace_back(AtomicNumber(6 ), 36.67100000,  34.29600000, 27.08900000);
  mol.emplace_back(AtomicNumber(8 ), 37.30500000,  33.23300000, 26.79500000);
  mol.emplace_back(AtomicNumber(6 ), 34.71700000,  32.76000000, 27.28600000);
  mol.emplace_back(AtomicNumber(6 ), 35.75200000,  32.05400000, 28.16000000);
  mol.emplace_back(AtomicNumber(6 ), 35.61200000,  30.57700000, 28.04400000);
  mol.emplace_back(AtomicNumber(7 ), 35.04000000,  30.25200000, 26.73000000);
  mol.emplace_back(AtomicNumber(6 ), 34.33800000,  29.10300000, 26.65000000);
  mol.emplace_back(AtomicNumber(7 ), 34.11000000,  28.43700000, 27.76800000);
  mol.emplace_back(AtomicNumber(7 ), 34.01400000,  28.65700000, 25.45700000);
  mol.emplace_back(AtomicNumber(7 ), 37.19700000,  35.39700000, 27.51300000);
  mol.emplace_back(AtomicNumber(6 ), 38.66800000,  35.50200000, 27.68000000);
  mol.emplace_back(AtomicNumber(6 ), 39.07600000,  34.93100000, 29.03100000);
  mol.emplace_back(AtomicNumber(8 ), 38.29700000,  34.94600000, 29.99600000);
  mol.emplace_back(AtomicNumber(6 ), 39.08000000,  36.94100000, 27.40600000);
  mol.emplace_back(AtomicNumber(6 ), 39.50200000,  37.34000000, 26.00200000);
  mol.emplace_back(AtomicNumber(6 ), 38.68400000,  36.64700000, 24.92300000);
  mol.emplace_back(AtomicNumber(6 ), 39.33700000,  38.85400000, 25.86200000);
  mol.emplace_back(AtomicNumber(7 ), 40.29400000,  34.41200000, 29.04500000);
  mol.emplace_back(AtomicNumber(6 ), 40.87300000,  33.80200000, 30.25300000);
  mol.emplace_back(AtomicNumber(6 ), 41.76500000,  34.82900000, 30.94400000);
  mol.emplace_back(AtomicNumber(8 ), 42.94500000,  34.99400000, 30.58300000);
  mol.emplace_back(AtomicNumber(6 ), 41.65100000,  32.52900000, 29.92300000);
  mol.emplace_back(AtomicNumber(6 ), 41.60800000,  31.44400000, 30.98900000);
  mol.emplace_back(AtomicNumber(6 ), 41.89600000,  30.08000000, 30.45600000);
  mol.emplace_back(AtomicNumber(7 ), 43.31100000,  29.73500000, 30.56300000);
  mol.emplace_back(AtomicNumber(6 ), 44.17400000,  29.90500000, 29.55400000);
  mol.emplace_back(AtomicNumber(7 ), 43.75400000,  30.31200000, 28.35600000);
  mol.emplace_back(AtomicNumber(7 ), 45.47700000,  29.72600000, 29.76300000);
  mol.emplace_back(AtomicNumber(7 ), 41.16500000,  35.53100000, 31.89800000);
  mol.emplace_back(AtomicNumber(6 ), 41.84500000,  36.55000000, 32.68600000);
  mol.emplace_back(AtomicNumber(6 ), 41.25100000,  37.94100000, 32.58800000);
  mol.emplace_back(AtomicNumber(8 ), 41.10200000,  38.52300000, 31.50000000);
  mol.emplace_back(AtomicNumber(7 ), 40.94600000,  38.47200000, 33.75700000);
  mol.emplace_back(AtomicNumber(6 ), 40.37300000,  39.81300000, 33.94400000);
  mol.emplace_back(AtomicNumber(6 ), 40.03100000,  39.99200000, 35.43200000);
  mol.emplace_back(AtomicNumber(8 ), 38.93300000,  40.52500000, 35.68700000);
  mol.emplace_back(AtomicNumber(8 ), 40.86200000,  39.57500000, 36.25100000);
  mol.emplace_back(AtomicNumber(1 ), 26.92359000,  23.39405000,  2.70240000);
  mol.emplace_back(AtomicNumber(1 ), 27.76664000,  24.57449000,  1.58857000);
  mol.emplace_back(AtomicNumber(1 ), 25.80382000,  25.69137000,  1.83776000);
  mol.emplace_back(AtomicNumber(1 ), 24.20066000,  25.53341000,  3.44380000);
  mol.emplace_back(AtomicNumber(1 ), 24.96122000,  23.79665000,  3.32771000);
  mol.emplace_back(AtomicNumber(1 ), 26.32946000,  24.31871000,  5.36451000);
  mol.emplace_back(AtomicNumber(1 ), 25.43351000,  25.91971000,  5.54648000);
  mol.emplace_back(AtomicNumber(1 ), 24.56109000,  22.91480000,  7.99868000);
  mol.emplace_back(AtomicNumber(1 ), 25.44564000,  24.52674000,  7.70799000);
  mol.emplace_back(AtomicNumber(1 ), 23.66291000,  24.52649000,  8.24491000);
  mol.emplace_back(AtomicNumber(1 ), 25.49362000,  27.80698000,  2.57791000);
  mol.emplace_back(AtomicNumber(1 ), 27.95299000,  28.94312000,  4.17541000);
  mol.emplace_back(AtomicNumber(1 ), 27.53638000,  29.98611000,  2.11255000);
  mol.emplace_back(AtomicNumber(1 ), 25.65555000,  30.09183000,  2.53684000);
  mol.emplace_back(AtomicNumber(1 ), 26.05602000,  31.75768000,  4.16567000);
  mol.emplace_back(AtomicNumber(1 ), 27.90909000,  31.64558000,  3.89354000);
  mol.emplace_back(AtomicNumber(1 ), 24.73835000,  32.17555000,  2.23338000);
  mol.emplace_back(AtomicNumber(1 ), 25.38196000,  33.43067000,  0.99814000);
  mol.emplace_back(AtomicNumber(1 ), 27.92525000,  29.69093000,  6.10622000);
  mol.emplace_back(AtomicNumber(1 ), 25.10614000,  30.12242000,  7.35162000);
  mol.emplace_back(AtomicNumber(1 ), 25.90724000,  29.49559000,  9.59909000);
  mol.emplace_back(AtomicNumber(1 ), 28.23549000,  28.03714000,  8.21591000);
  mol.emplace_back(AtomicNumber(1 ), 28.37174000,  29.73983000,  9.01684000);
  mol.emplace_back(AtomicNumber(1 ), 24.38117000,  28.02188000,  8.35722000);
  mol.emplace_back(AtomicNumber(1 ), 25.73798000,  27.43439000,  7.22620000);
  mol.emplace_back(AtomicNumber(1 ), 25.73786000,  26.93282000,  9.01920000);
  mol.emplace_back(AtomicNumber(1 ), 26.93576000,  27.88169000, 10.85744000);
  mol.emplace_back(AtomicNumber(1 ), 28.53974000,  27.10634000, 10.31767000);
  mol.emplace_back(AtomicNumber(1 ), 28.53947000,  28.79390000, 11.10415000);
  mol.emplace_back(AtomicNumber(1 ), 25.29979000,  31.70614000,  9.19908000);
  mol.emplace_back(AtomicNumber(1 ), 27.70543000,  33.71191000,  8.60355000);
  mol.emplace_back(AtomicNumber(1 ), 24.77323000,  34.22993000,  9.56088000);
  mol.emplace_back(AtomicNumber(1 ), 26.13797000,  35.50389000,  9.24860000);
  mol.emplace_back(AtomicNumber(1 ), 23.50686000,  33.40666000,  7.73847000);
  mol.emplace_back(AtomicNumber(1 ), 27.02383000,  35.86827000,  7.03026000);
  mol.emplace_back(AtomicNumber(1 ), 22.93086000,  33.50109000,  5.28193000);
  mol.emplace_back(AtomicNumber(1 ), 26.47070000,  35.89506000,  4.54332000);
  mol.emplace_back(AtomicNumber(1 ), 24.31321000,  34.90166000,  3.80234000);
  mol.emplace_back(AtomicNumber(1 ), 28.92666000,  34.41475000, 10.38549000);
  mol.emplace_back(AtomicNumber(1 ), 27.83983000,  33.35771000, 13.09064000);
  mol.emplace_back(AtomicNumber(1 ), 30.79188000,  33.86886000, 12.25902000);
  mol.emplace_back(AtomicNumber(1 ), 30.28238000,  32.34016000, 14.73286000);
  mol.emplace_back(AtomicNumber(1 ), 31.18697000,  33.95014000, 14.49841000);
  mol.emplace_back(AtomicNumber(1 ), 29.35037000,  33.94988000, 14.80393000);
  mol.emplace_back(AtomicNumber(1 ), 28.88092000,  31.47874000, 12.23978000);
  mol.emplace_back(AtomicNumber(1 ), 30.49302000,  31.75096000, 11.34965000);
  mol.emplace_back(AtomicNumber(1 ), 30.49277000,  31.19600000, 13.12685000);
  mol.emplace_back(AtomicNumber(1 ), 27.04413000,  35.16220000, 14.15881000);
  mol.emplace_back(AtomicNumber(1 ), 28.15111000,  38.00941000, 13.36474000);
  mol.emplace_back(AtomicNumber(1 ), 25.68905000,  37.41077000, 13.33535000);
  mol.emplace_back(AtomicNumber(1 ), 25.87828000,  37.10794000, 15.22986000);
  mol.emplace_back(AtomicNumber(1 ), 26.31729000,  39.42637000, 15.62953000);
  mol.emplace_back(AtomicNumber(1 ), 26.33094000,  39.81892000, 13.81653000);
  mol.emplace_back(AtomicNumber(1 ), 23.92768000,  38.84916000, 13.68113000);
  mol.emplace_back(AtomicNumber(1 ), 23.98424000,  38.75992000, 15.58003000);
  mol.emplace_back(AtomicNumber(1 ), 24.34491000,  41.22267000, 15.65433000);
  mol.emplace_back(AtomicNumber(1 ), 24.15276000,  41.29991000, 13.80881000);
  mol.emplace_back(AtomicNumber(1 ), 22.00532000,  39.66350000, 14.94620000);
  mol.emplace_back(AtomicNumber(1 ), 22.09926000,  41.24732000, 15.85584000);
  mol.emplace_back(AtomicNumber(1 ), 29.60442000,  39.03937000, 14.56928000);
  mol.emplace_back(AtomicNumber(1 ), 30.16374000,  37.71144000, 17.31625000);
  mol.emplace_back(AtomicNumber(1 ), 32.42146000,  38.65951000, 17.18921000);
  mol.emplace_back(AtomicNumber(1 ), 30.69834000,  40.59943000, 15.70596000);
  mol.emplace_back(AtomicNumber(1 ), 31.37527000,  37.51759000, 14.61145000);
  mol.emplace_back(AtomicNumber(1 ), 32.89377000,  38.58195000, 14.44835000);
  mol.emplace_back(AtomicNumber(1 ), 32.89347000,  37.13293000, 15.61743000);
  mol.emplace_back(AtomicNumber(1 ), 30.91428000,  39.43490000, 19.03919000);
  mol.emplace_back(AtomicNumber(1 ), 28.47380000,  41.23734000, 19.35672000);
  mol.emplace_back(AtomicNumber(1 ), 31.03644000,  40.67897000, 21.01796000);
  mol.emplace_back(AtomicNumber(1 ), 29.56406000,  41.81326000, 21.50476000);
  mol.emplace_back(AtomicNumber(1 ), 29.67267000,  38.77855000, 21.16901000);
  mol.emplace_back(AtomicNumber(1 ), 29.34094000,  38.68258000, 23.49392000);
  mol.emplace_back(AtomicNumber(1 ), 30.26609000,  40.29441000, 23.38705000);
  mol.emplace_back(AtomicNumber(1 ), 28.41055000,  40.29415000, 23.54007000);
  mol.emplace_back(AtomicNumber(1 ), 27.33049000,  38.64614000, 21.14004000);
  mol.emplace_back(AtomicNumber(1 ), 27.11458000,  40.25825000, 22.04546000);
  mol.emplace_back(AtomicNumber(1 ), 27.54588000,  40.25799000, 20.23428000);
  mol.emplace_back(AtomicNumber(1 ), 31.39159000,  41.66869000, 17.55397000);
  mol.emplace_back(AtomicNumber(1 ), 31.10054000,  44.72865000, 18.31082000);
  mol.emplace_back(AtomicNumber(1 ), 33.43540000,  44.95987000, 17.36873000);
  mol.emplace_back(AtomicNumber(1 ), 32.19506000,  42.77355000, 15.85810000);
  mol.emplace_back(AtomicNumber(1 ), 33.62203000,  42.18130000, 18.60871000);
  mol.emplace_back(AtomicNumber(1 ), 34.74279000,  43.66401000, 18.71242000);
  mol.emplace_back(AtomicNumber(1 ), 33.11306000,  43.66385000, 19.61265000);
  mol.emplace_back(AtomicNumber(1 ), 29.39472000,  42.84264000, 16.49806000);
  mol.emplace_back(AtomicNumber(1 ), 27.96053000,  43.45359000, 14.76695000);
  mol.emplace_back(AtomicNumber(1 ), 28.93135000,  45.09707000, 14.61098000);
  mol.emplace_back(AtomicNumber(1 ), 30.90273000,  42.33463000, 14.48152000);
  mol.emplace_back(AtomicNumber(1 ), 31.11309000,  42.73542000, 11.45339000);
  mol.emplace_back(AtomicNumber(1 ), 33.21916000,  42.69063000, 12.73358000);
  mol.emplace_back(AtomicNumber(1 ), 32.79019000,  40.97394000, 13.36144000);
  mol.emplace_back(AtomicNumber(1 ), 32.73691000,  40.11319000, 10.98556000);
  mol.emplace_back(AtomicNumber(1 ), 33.16585000,  41.80581000, 10.35039000);
  mol.emplace_back(AtomicNumber(1 ), 35.15807000,  41.57258000, 12.20777000);
  mol.emplace_back(AtomicNumber(1 ), 34.82463000,  39.73184000, 11.86604000);
  mol.emplace_back(AtomicNumber(1 ), 35.70645000,  41.93590000,  9.91543000);
  mol.emplace_back(AtomicNumber(1 ), 36.61202000,  40.37457000, 10.52351000);
  mol.emplace_back(AtomicNumber(1 ), 33.98016000,  40.08364000,  9.11665000);
  mol.emplace_back(AtomicNumber(1 ), 35.46423000,  40.52906000,  8.14510000);
  mol.emplace_back(AtomicNumber(1 ), 30.36882000,  41.01286000, 10.06502000);
  mol.emplace_back(AtomicNumber(1 ), 29.38652000,  38.52629000, 11.66872000);
  mol.emplace_back(AtomicNumber(1 ), 28.11525000,  39.58633000,  9.00958000);
  mol.emplace_back(AtomicNumber(1 ), 26.75838000,  38.97330000, 11.61030000);
  mol.emplace_back(AtomicNumber(1 ), 26.44804000,  37.64398000,  9.71690000);
  mol.emplace_back(AtomicNumber(1 ), 27.96353000,  37.29942000,  8.69228000);
  mol.emplace_back(AtomicNumber(1 ), 27.96335000,  36.91325000, 10.51362000);
  mol.emplace_back(AtomicNumber(1 ), 30.39707000,  36.56232000, 10.99419000);
  mol.emplace_back(AtomicNumber(1 ), 32.11862000,  36.94713000,  8.33484000);
  mol.emplace_back(AtomicNumber(1 ), 33.79734000,  36.46047000,  9.36622000);
  mol.emplace_back(AtomicNumber(1 ), 32.40241000,  33.78680000, 10.20063000);
  mol.emplace_back(AtomicNumber(1 ), 33.53938000,  34.11832000,  8.75550000);
  mol.emplace_back(AtomicNumber(1 ), 32.06453000,  36.52485000, 11.86858000);
  mol.emplace_back(AtomicNumber(1 ), 33.67242000,  35.60469000, 12.05104000);
  mol.emplace_back(AtomicNumber(1 ), 33.67215000,  37.37175000, 11.46459000);
  mol.emplace_back(AtomicNumber(1 ), 34.63760000,  32.89501000, 10.79346000);
  mol.emplace_back(AtomicNumber(1 ), 35.47704000,  34.46503000, 10.24971000);
  mol.emplace_back(AtomicNumber(1 ), 34.39293000,  34.46482000, 11.76336000);
  mol.emplace_back(AtomicNumber(1 ), 31.98727000,  35.57532000,  6.67540000);
  mol.emplace_back(AtomicNumber(1 ), 29.55891000,  33.58632000,  7.07405000);
  mol.emplace_back(AtomicNumber(1 ), 29.34435000,  33.72637000,  4.53236000);
  mol.emplace_back(AtomicNumber(1 ), 31.85235000,  33.50934000,  4.30635000);
  mol.emplace_back(AtomicNumber(1 ), 28.28779000,  35.62304000,  5.11500000);
  mol.emplace_back(AtomicNumber(1 ), 29.73455000,  36.30590000,  6.06677000);
  mol.emplace_back(AtomicNumber(1 ), 29.73441000,  36.36161000,  4.20577000);
  mol.emplace_back(AtomicNumber(1 ), 29.82418000,  31.41032000,  6.87340000);
  mol.emplace_back(AtomicNumber(1 ), 32.75565000,  30.47283000,  6.32757000);
  mol.emplace_back(AtomicNumber(1 ), 30.81835000,  30.36759000,  8.57605000);
  mol.emplace_back(AtomicNumber(1 ), 31.37899000,  28.60098000,  7.74696000);
  mol.emplace_back(AtomicNumber(1 ), 33.11442000,  28.44188000,  8.81167000);
  mol.emplace_back(AtomicNumber(1 ), 34.78193000,  29.78824000,  9.00545000);
  mol.emplace_back(AtomicNumber(1 ), 33.75151000,  31.04703000,  8.10048000);
  mol.emplace_back(AtomicNumber(1 ), 33.79873000,  31.04672000,  9.96171000);
  mol.emplace_back(AtomicNumber(1 ), 31.77366000,  28.37563000, 10.77969000);
  mol.emplace_back(AtomicNumber(1 ), 32.60932000,  29.98773000, 11.18969000);
  mol.emplace_back(AtomicNumber(1 ), 30.94316000,  29.98747000, 10.35883000);
  mol.emplace_back(AtomicNumber(1 ), 32.90335000,  28.46009000,  5.28622000);
  mol.emplace_back(AtomicNumber(1 ), 30.11772000,  27.53898000,  4.06195000);
  mol.emplace_back(AtomicNumber(1 ), 32.96170000,  27.31028000,  2.99262000);
  mol.emplace_back(AtomicNumber(1 ), 31.48188000,  26.28174000,  2.42540000);
  mol.emplace_back(AtomicNumber(1 ), 30.23091000,  28.49528000,  2.06538000);
  mol.emplace_back(AtomicNumber(1 ), 31.95851000,  29.31860000,  2.30757000);
  mol.emplace_back(AtomicNumber(1 ), 33.47809000,  27.62691000,  1.09204000);
  mol.emplace_back(AtomicNumber(1 ), 29.37434000,  25.86344000,  5.02073000);
  mol.emplace_back(AtomicNumber(1 ), 31.33133000,  23.78687000,  6.15857000);
  mol.emplace_back(AtomicNumber(1 ), 29.96919000,  23.61043000,  8.24436000);
  mol.emplace_back(AtomicNumber(1 ), 31.62084000,  24.86706000,  8.96189000);
  mol.emplace_back(AtomicNumber(1 ), 31.40665000,  26.15387000,  7.63390000);
  mol.emplace_back(AtomicNumber(1 ), 30.28485000,  26.15356000,  9.11984000);
  mol.emplace_back(AtomicNumber(1 ), 27.73750000,  24.26679000,  7.72570000);
  mol.emplace_back(AtomicNumber(1 ), 28.33897000,  25.77263000,  8.64001000);
  mol.emplace_back(AtomicNumber(1 ), 28.28688000,  25.77233000,  6.77891000);
  mol.emplace_back(AtomicNumber(1 ), 30.18195000,  21.91157000,  6.94441000);
  mol.emplace_back(AtomicNumber(1 ), 27.71941000,  21.23469000,  5.17231000);
  mol.emplace_back(AtomicNumber(1 ), 30.10868000,  19.53909000,  6.19333000);
  mol.emplace_back(AtomicNumber(1 ), 28.45778000,  18.84368000,  5.53892000);
  mol.emplace_back(AtomicNumber(1 ), 30.53820000,  20.55607000,  4.02176000);
  mol.emplace_back(AtomicNumber(1 ), 30.06964000,  18.69004000,  3.83922000);
  mol.emplace_back(AtomicNumber(1 ), 30.39755000,  20.47498000,  1.86597000);
  mol.emplace_back(AtomicNumber(1 ), 25.72896000,  20.68215000,  9.23891000);
  mol.emplace_back(AtomicNumber(1 ), 24.49907000,  18.21589000,  7.72180000);
  mol.emplace_back(AtomicNumber(1 ), 23.61320000,  19.46455000,  8.69337000);
  mol.emplace_back(AtomicNumber(1 ), 23.53275000,  19.65343000,  6.01009000);
  mol.emplace_back(AtomicNumber(1 ), 23.91498000,  21.18953000,  6.89890000);
  mol.emplace_back(AtomicNumber(1 ), 25.87519000,  20.98813000,  5.36661000);
  mol.emplace_back(AtomicNumber(1 ), 25.92870000,  19.19139000,  5.43009000);
  mol.emplace_back(AtomicNumber(1 ), 27.52378000,  18.15403000,  7.50673000);
  mol.emplace_back(AtomicNumber(1 ), 27.30865000,  16.34801000,  9.92196000);
  mol.emplace_back(AtomicNumber(1 ), 29.39732000,  15.23212000,  8.61738000);
  mol.emplace_back(AtomicNumber(1 ), 27.64570000,  15.13396000,  7.93848000);
  mol.emplace_back(AtomicNumber(1 ), 28.06407000,  16.74572000,  6.38239000);
  mol.emplace_back(AtomicNumber(1 ), 29.02674000,  19.25967000,  9.18937000);
  mol.emplace_back(AtomicNumber(1 ), 31.66541000,  18.36509000, 10.39767000);
  mol.emplace_back(AtomicNumber(1 ), 30.19349000,  21.10307000,  9.87689000);
  mol.emplace_back(AtomicNumber(1 ), 31.80316000,  21.03532000, 10.82823000);
  mol.emplace_back(AtomicNumber(1 ), 32.43463000,  18.53692000,  9.10939000);
  mol.emplace_back(AtomicNumber(1 ), 32.46185000,  18.64472000, 12.42663000);
  mol.emplace_back(AtomicNumber(1 ), 30.37318000,  18.70166000, 14.62958000);
  mol.emplace_back(AtomicNumber(1 ), 32.39984000,  18.24026000, 16.13029000);
  mol.emplace_back(AtomicNumber(1 ), 34.45693000,  17.91129000, 14.26312000);
  mol.emplace_back(AtomicNumber(1 ), 31.24730000,  16.42656000, 14.56315000);
  mol.emplace_back(AtomicNumber(1 ), 32.85282000,  16.50621000, 13.62440000);
  mol.emplace_back(AtomicNumber(1 ), 32.85258000,  16.06426000, 15.43302000);
  mol.emplace_back(AtomicNumber(1 ), 30.58583000,  20.11137000, 16.43355000);
  mol.emplace_back(AtomicNumber(1 ), 30.85554000,  22.97929000, 15.70507000);
  mol.emplace_back(AtomicNumber(1 ), 30.82691000,  21.53235000, 18.54894000);
  mol.emplace_back(AtomicNumber(1 ), 28.89776000,  20.96265000, 17.09454000);
  mol.emplace_back(AtomicNumber(1 ), 28.43508000,  22.27196000, 18.40181000);
  mol.emplace_back(AtomicNumber(1 ), 29.85586000,  24.24102000, 18.58814000);
  mol.emplace_back(AtomicNumber(1 ), 31.43490000,  23.57825000, 19.31796000);
  mol.emplace_back(AtomicNumber(1 ), 31.43462000,  24.30577000, 17.60416000);
  mol.emplace_back(AtomicNumber(1 ), 27.26922000,  22.90224000, 16.42041000);
  mol.emplace_back(AtomicNumber(1 ), 28.75260000,  24.02089000, 16.53646000);
  mol.emplace_back(AtomicNumber(1 ), 28.75244000,  22.59587000, 15.33825000);
  mol.emplace_back(AtomicNumber(1 ), 33.11293000,  20.56021000, 17.17309000);
  mol.emplace_back(AtomicNumber(1 ), 35.24144000,  22.51151000, 17.86398000);
  mol.emplace_back(AtomicNumber(1 ), 35.05285000,  19.96691000, 18.31261000);
  mol.emplace_back(AtomicNumber(1 ), 35.66193000,  19.81815000, 16.45679000);
  mol.emplace_back(AtomicNumber(1 ), 37.77881000,  20.76053000, 17.06116000);
  mol.emplace_back(AtomicNumber(1 ), 37.24723000,  20.84900000, 18.86005000);
  mol.emplace_back(AtomicNumber(1 ), 37.39507000,  19.35485000, 20.02396000);
  mol.emplace_back(AtomicNumber(1 ), 34.36591000,  20.87468000, 14.77448000);
  mol.emplace_back(AtomicNumber(1 ), 36.72561000,  21.84546000, 13.31092000);
  mol.emplace_back(AtomicNumber(1 ), 35.46760000,  19.91345000, 12.47342000);
  mol.emplace_back(AtomicNumber(1 ), 33.92490000,  20.92354000, 12.28555000);
  mol.emplace_back(AtomicNumber(1 ), 33.58221000,  21.41441000, 10.13335000);
  mol.emplace_back(AtomicNumber(1 ), 34.93077000,  21.91964000,  8.93271000);
  mol.emplace_back(AtomicNumber(1 ), 33.34699000,  23.03113000, 13.72601000);
  mol.emplace_back(AtomicNumber(1 ), 33.58488000,  25.30295000, 11.85796000);
  mol.emplace_back(AtomicNumber(1 ), 31.91004000,  24.86424000, 14.42365000);
  mol.emplace_back(AtomicNumber(1 ), 30.54898000,  26.85542000, 13.34043000);
  mol.emplace_back(AtomicNumber(1 ), 32.11239000,  27.26235000, 14.26534000);
  mol.emplace_back(AtomicNumber(1 ), 32.11218000,  27.23467000, 12.40371000);
  mol.emplace_back(AtomicNumber(1 ), 31.34748000,  23.31424000, 12.66795000);
  mol.emplace_back(AtomicNumber(1 ), 30.06315000,  24.66095000, 12.71456000);
  mol.emplace_back(AtomicNumber(1 ), 31.36438000,  24.66084000, 11.38294000);
  mol.emplace_back(AtomicNumber(1 ), 34.33276000,  24.92679000, 15.40663000);
  mol.emplace_back(AtomicNumber(1 ), 35.12085000,  27.75118000, 15.74831000);
  mol.emplace_back(AtomicNumber(1 ), 35.95607000,  25.08880000, 17.16512000);
  mol.emplace_back(AtomicNumber(1 ), 36.55033000,  26.75467000, 17.71742000);
  mol.emplace_back(AtomicNumber(1 ), 34.04674000,  27.54210000, 17.80746000);
  mol.emplace_back(AtomicNumber(1 ), 33.52741000,  25.80677000, 17.42147000);
  mol.emplace_back(AtomicNumber(1 ), 34.79771000,  24.97703000, 19.43951000);
  mol.emplace_back(AtomicNumber(1 ), 35.32444000,  26.73128000, 19.81452000);
  mol.emplace_back(AtomicNumber(1 ), 32.91398000,  27.41174000, 20.07018000);
  mol.emplace_back(AtomicNumber(1 ), 32.35446000,  25.70256000, 19.67003000);
  mol.emplace_back(AtomicNumber(1 ), 32.47594000,  25.76257000, 22.06955000);
  mol.emplace_back(AtomicNumber(1 ), 34.03769000,  24.95376000, 21.56800000);
  mol.emplace_back(AtomicNumber(1 ), 36.95680000,  24.80880000, 14.64239000);
  mol.emplace_back(AtomicNumber(1 ), 39.54641000,  26.24584000, 14.58598000);
  mol.emplace_back(AtomicNumber(1 ), 38.38786000,  23.64910000, 13.41464000);
  mol.emplace_back(AtomicNumber(1 ), 39.91558000,  24.35211000, 12.61641000);
  mol.emplace_back(AtomicNumber(1 ), 39.91528000,  23.95155000, 14.43465000);
  mol.emplace_back(AtomicNumber(1 ), 36.83623000,  25.86488000, 12.14471000);
  mol.emplace_back(AtomicNumber(1 ), 38.36655000,  27.18126000,  9.99450000);
  mol.emplace_back(AtomicNumber(1 ), 35.30541000,  27.22417000, 10.60683000);
  mol.emplace_back(AtomicNumber(1 ), 36.15223000,  27.73828000,  8.99713000);
  mol.emplace_back(AtomicNumber(1 ), 36.90418000,  25.49616000,  8.56048000);
  mol.emplace_back(AtomicNumber(1 ), 36.41732000,  24.91834000, 10.26772000);
  mol.emplace_back(AtomicNumber(1 ), 34.14751000,  24.84931000,  9.76105000);
  mol.emplace_back(AtomicNumber(1 ), 34.28825000,  26.23774000,  8.49791000);
  mol.emplace_back(AtomicNumber(1 ), 35.68707000,  23.57220000,  7.91414000);
  mol.emplace_back(AtomicNumber(1 ), 33.82073000,  23.66898000,  7.78262000);
  mol.emplace_back(AtomicNumber(1 ), 33.89196000,  25.15303000,  6.03339000);
  mol.emplace_back(AtomicNumber(1 ), 35.47054000,  25.90795000,  6.56533000);
  mol.emplace_back(AtomicNumber(1 ), 36.36146000,  28.38045000, 12.78050000);
  mol.emplace_back(AtomicNumber(1 ), 36.34707000,  31.22614000, 11.79548000);
  mol.emplace_back(AtomicNumber(1 ), 35.93784000,  29.94964000, 14.55694000);
  mol.emplace_back(AtomicNumber(1 ), 33.93637000,  31.54836000, 12.70976000);
  mol.emplace_back(AtomicNumber(1 ), 34.13217000,  29.69703000, 12.67094000);
  mol.emplace_back(AtomicNumber(1 ), 34.84656000,  32.57652000, 14.73931000);
  mol.emplace_back(AtomicNumber(1 ), 36.44924000,  31.99360000, 15.48557000);
  mol.emplace_back(AtomicNumber(1 ), 36.44897000,  32.85028000, 13.83254000);
  mol.emplace_back(AtomicNumber(1 ), 33.10626000,  29.40213000, 14.78416000);
  mol.emplace_back(AtomicNumber(1 ), 33.76648000,  31.00901000, 15.45298000);
  mol.emplace_back(AtomicNumber(1 ), 32.28829000,  31.00878000, 14.32101000);
  mol.emplace_back(AtomicNumber(1 ), 38.49805000,  29.12618000, 13.94904000);
  mol.emplace_back(AtomicNumber(1 ), 40.17054000,  31.41708000, 14.79579000);
  mol.emplace_back(AtomicNumber(1 ), 40.44643000,  29.17620000, 15.83862000);
  mol.emplace_back(AtomicNumber(1 ), 41.02373000,  28.46041000, 14.21857000);
  mol.emplace_back(AtomicNumber(1 ), 43.04148000,  29.92492000, 14.21355000);
  mol.emplace_back(AtomicNumber(1 ), 42.42857000,  30.86636000, 15.69413000);
  mol.emplace_back(AtomicNumber(1 ), 44.08772000,  30.31121000, 17.16781000);
  mol.emplace_back(AtomicNumber(1 ), 44.37047000,  28.54027000, 17.71485000);
  mol.emplace_back(AtomicNumber(1 ), 40.39728000,  28.99838000, 12.11427000);
  mol.emplace_back(AtomicNumber(1 ), 42.84673000,  30.07171000, 10.79498000);
  mol.emplace_back(AtomicNumber(1 ), 41.22400000,  27.87941000, 10.48695000);
  mol.emplace_back(AtomicNumber(1 ), 40.45020000,  28.91641000,  9.19140000);
  mol.emplace_back(AtomicNumber(1 ), 43.65514000,  28.18170000, 10.59941000);
  mol.emplace_back(AtomicNumber(1 ), 39.35191000,  31.20482000, 10.52579000);
  mol.emplace_back(AtomicNumber(1 ), 40.71240000,  33.16357000,  8.56001000);
  mol.emplace_back(AtomicNumber(1 ), 38.68005000,  33.62344000,  7.53326000);
  mol.emplace_back(AtomicNumber(1 ), 38.73476000,  31.71351000,  7.99966000);
  mol.emplace_back(AtomicNumber(1 ), 36.82738000,  31.75539000,  8.91532000);
  mol.emplace_back(AtomicNumber(1 ), 37.18654000,  33.31411000,  9.85477000);
  mol.emplace_back(AtomicNumber(1 ), 35.28776000,  33.67531000,  8.24451000);
  mol.emplace_back(AtomicNumber(1 ), 36.79727000,  34.65521000,  7.71864000);
  mol.emplace_back(AtomicNumber(1 ), 37.37686000,  32.63851000,  6.17116000);
  mol.emplace_back(AtomicNumber(1 ), 35.67353000,  32.00681000,  6.50460000);
  mol.emplace_back(AtomicNumber(1 ), 34.64868000,  33.91719000,  5.51645000);
  mol.emplace_back(AtomicNumber(1 ), 36.12255000,  34.97226000,  5.75997000);
  mol.emplace_back(AtomicNumber(1 ), 39.47354000,  33.40504000, 11.80850000);
  mol.emplace_back(AtomicNumber(1 ), 39.97639000,  36.36539000, 11.33742000);
  mol.emplace_back(AtomicNumber(1 ), 38.03160000,  34.93533000, 13.37680000);
  mol.emplace_back(AtomicNumber(1 ), 38.37457000,  36.80332000, 13.25809000);
  mol.emplace_back(AtomicNumber(1 ), 37.25691000,  35.18752000, 10.87965000);
  mol.emplace_back(AtomicNumber(1 ), 36.16539000,  35.87921000, 12.24216000);
  mol.emplace_back(AtomicNumber(1 ), 36.59811000,  36.39212000,  9.37372000);
  mol.emplace_back(AtomicNumber(1 ), 41.05746000,  33.52137000, 12.84488000);
  mol.emplace_back(AtomicNumber(1 ), 42.98955000,  33.37521000, 14.13894000);
  mol.emplace_back(AtomicNumber(1 ), 42.85755000,  35.28726000, 14.42603000);
  mol.emplace_back(AtomicNumber(1 ), 39.90869000,  33.93397000, 15.31592000);
  mol.emplace_back(AtomicNumber(1 ), 40.84630000,  34.36665000, 18.21005000);
  mol.emplace_back(AtomicNumber(1 ), 38.09788000,  33.39065000, 16.93644000);
  mol.emplace_back(AtomicNumber(1 ), 38.77262000,  36.30571000, 17.83964000);
  mol.emplace_back(AtomicNumber(1 ), 39.10735000,  35.69538000, 16.11100000);
  mol.emplace_back(AtomicNumber(1 ), 38.05432000,  32.82133000, 19.26700000);
  mol.emplace_back(AtomicNumber(1 ), 38.85921000,  34.42322000, 19.76837000);
  mol.emplace_back(AtomicNumber(1 ), 37.12002000,  34.42299000, 19.10381000);
  mol.emplace_back(AtomicNumber(1 ), 36.45795000,  34.72865000, 16.55937000);
  mol.emplace_back(AtomicNumber(1 ), 36.44005000,  36.32753000, 17.51256000);
  mol.emplace_back(AtomicNumber(1 ), 36.87682000,  36.32725000, 15.70268000);
  mol.emplace_back(AtomicNumber(1 ), 42.04062000,  30.05628000, 18.90046000);
  mol.emplace_back(AtomicNumber(1 ), 41.45001000,  31.31552000, 21.75827000);
  mol.emplace_back(AtomicNumber(1 ), 42.90170000,  30.32921000, 21.28655000);
  mol.emplace_back(AtomicNumber(1 ), 43.24801000,  32.97326000, 21.48492000);
  mol.emplace_back(AtomicNumber(1 ), 43.85093000,  32.25174000, 19.93933000);
  mol.emplace_back(AtomicNumber(1 ), 42.26648000,  33.98644000, 19.13143000);
  mol.emplace_back(AtomicNumber(1 ), 41.10605000,  33.76663000, 20.48547000);
  mol.emplace_back(AtomicNumber(1 ), 37.90072000,  28.49259000, 19.40249000);
  mol.emplace_back(AtomicNumber(1 ), 39.54776000,  25.94844000, 20.11677000);
  mol.emplace_back(AtomicNumber(1 ), 38.14695000,  26.09403000, 18.97278000);
  mol.emplace_back(AtomicNumber(1 ), 40.63029000,  25.95549000, 17.88287000);
  mol.emplace_back(AtomicNumber(1 ), 39.60461000,  27.28825000, 17.23673000);
  mol.emplace_back(AtomicNumber(1 ), 41.62800000,  28.51118000, 18.11393000);
  mol.emplace_back(AtomicNumber(1 ), 41.83459000,  27.43335000, 19.53780000);
  mol.emplace_back(AtomicNumber(1 ), 40.40606000,  28.14533000, 21.91815000);
  mol.emplace_back(AtomicNumber(1 ), 38.29365000,  27.25265000, 23.92090000);
  mol.emplace_back(AtomicNumber(1 ), 40.25444000,  26.61153000, 24.82136000);
  mol.emplace_back(AtomicNumber(1 ), 41.22700000,  27.82960000, 23.76523000);
  mol.emplace_back(AtomicNumber(1 ), 40.69219000,  30.17813000, 24.45890000);
  mol.emplace_back(AtomicNumber(1 ), 38.95785000,  30.30016000, 22.40498000);
  mol.emplace_back(AtomicNumber(1 ), 37.67056000,  31.63283000, 24.85000000);
  mol.emplace_back(AtomicNumber(1 ), 38.63765000,  32.77845000, 22.04980000);
  mol.emplace_back(AtomicNumber(1 ), 37.97577000,  33.81913000, 23.43334000);
  mol.emplace_back(AtomicNumber(1 ), 39.79055000,  32.90013000, 24.97180000);
  mol.emplace_back(AtomicNumber(1 ), 40.45747000,  31.90013000, 23.52411000);
  mol.emplace_back(AtomicNumber(1 ), 39.07331000,  34.94439000, 22.60165000);
  mol.emplace_back(AtomicNumber(1 ), 40.69351000,  35.87332000, 22.43600000);
  mol.emplace_back(AtomicNumber(1 ), 36.71024000,  30.07988000, 21.92956000);
  mol.emplace_back(AtomicNumber(1 ), 34.47852000,  31.98361000, 21.41596000);
  mol.emplace_back(AtomicNumber(1 ), 35.08765000,  29.12398000, 20.19247000);
  mol.emplace_back(AtomicNumber(1 ), 33.89858000,  30.42780000, 19.51115000);
  mol.emplace_back(AtomicNumber(1 ), 35.80962000,  31.97302000, 19.08021000);
  mol.emplace_back(AtomicNumber(1 ), 37.00372000,  30.71023000, 19.76290000);
  mol.emplace_back(AtomicNumber(1 ), 38.12144000,  30.54986000, 17.67336000);
  mol.emplace_back(AtomicNumber(1 ), 37.27750000,  29.64614000, 16.26393000);
  mol.emplace_back(AtomicNumber(1 ), 32.49930000,  31.95564000, 21.96472000);
  mol.emplace_back(AtomicNumber(1 ), 31.24195000,  29.21669000, 23.02621000);
  mol.emplace_back(AtomicNumber(1 ), 31.51492000,  30.43839000, 24.89955000);
  mol.emplace_back(AtomicNumber(1 ), 30.99908000,  32.05210000, 23.97798000);
  mol.emplace_back(AtomicNumber(1 ), 28.72920000,  31.39654000, 23.91044000);
  mol.emplace_back(AtomicNumber(1 ), 29.12112000,  29.69389000, 24.57943000);
  mol.emplace_back(AtomicNumber(1 ), 29.87237000,  30.64578000, 26.71702000);
  mol.emplace_back(AtomicNumber(1 ), 29.65953000,  32.38527000, 26.04816000);
  mol.emplace_back(AtomicNumber(1 ), 27.41139000,  30.36500000, 26.66787000);
  mol.emplace_back(AtomicNumber(1 ), 28.68296000,  33.79754000, 26.25354000);
  mol.emplace_back(AtomicNumber(1 ), 26.98098000,  34.51458000, 26.57680000);
  mol.emplace_back(AtomicNumber(1 ), 25.35718000,  31.31597000, 26.75236000);
  mol.emplace_back(AtomicNumber(1 ), 29.95783000,  28.44894000, 21.41165000);
  mol.emplace_back(AtomicNumber(1 ), 28.83107000,  30.64086000, 19.51293000);
  mol.emplace_back(AtomicNumber(1 ), 29.23828000,  27.60395000, 19.18774000);
  mol.emplace_back(AtomicNumber(1 ), 28.33072000,  28.83913000, 17.98504000);
  mol.emplace_back(AtomicNumber(1 ), 31.30146000,  29.08525000, 18.67682000);
  mol.emplace_back(AtomicNumber(1 ), 29.76796000,  27.20329000, 16.77731000);
  mol.emplace_back(AtomicNumber(1 ), 31.34358000,  28.01792000, 16.21244000);
  mol.emplace_back(AtomicNumber(1 ), 31.34330000,  26.94640000, 17.73502000);
  mol.emplace_back(AtomicNumber(1 ), 29.08821000,  30.35749000, 17.00061000);
  mol.emplace_back(AtomicNumber(1 ), 30.68505000,  31.06289000, 17.64701000);
  mol.emplace_back(AtomicNumber(1 ), 30.68482000,  30.01837000, 16.10578000);
  mol.emplace_back(AtomicNumber(1 ), 26.76269000,  31.06036000, 19.38845000);
  mol.emplace_back(AtomicNumber(1 ), 24.91612000,  29.20436000, 20.99534000);
  mol.emplace_back(AtomicNumber(1 ), 24.67941000,  32.34059000, 20.60667000);
  mol.emplace_back(AtomicNumber(1 ), 26.66599000,  31.89778000, 21.96567000);
  mol.emplace_back(AtomicNumber(1 ), 25.26852000,  32.54628000, 23.04490000);
  mol.emplace_back(AtomicNumber(1 ), 22.90068000,  30.20885000, 21.92744000);
  mol.emplace_back(AtomicNumber(1 ), 23.11000000,  31.81646000, 22.84235000);
  mol.emplace_back(AtomicNumber(1 ), 22.46565000,  31.81622000, 21.09557000);
  mol.emplace_back(AtomicNumber(1 ), 24.73960000,  30.06103000, 23.64012000);
  mol.emplace_back(AtomicNumber(1 ), 26.34582000,  29.60859000, 22.81511000);
  mol.emplace_back(AtomicNumber(1 ), 26.34555000,  30.73311000, 24.29898000);
  mol.emplace_back(AtomicNumber(1 ), 23.12305000,  28.49688000, 20.07480000);
  mol.emplace_back(AtomicNumber(1 ), 22.01681000,  30.05325000, 17.63064000);
  mol.emplace_back(AtomicNumber(1 ), 23.45886000,  28.40630000, 16.61238000);
  mol.emplace_back(AtomicNumber(1 ), 22.59647000,  27.05021000, 17.63385000);
  mol.emplace_back(AtomicNumber(1 ), 21.98691000,  29.69333000, 14.99136000);
  mol.emplace_back(AtomicNumber(1 ), 20.75313000,  25.99400000, 16.85982000);
  mol.emplace_back(AtomicNumber(1 ), 20.22395000,  29.41187000, 13.22077000);
  mol.emplace_back(AtomicNumber(1 ), 19.00295000,  25.66833000, 15.08083000);
  mol.emplace_back(AtomicNumber(1 ), 18.89726000,  27.27222000, 13.18909000);
  mol.emplace_back(AtomicNumber(1 ), 19.97402000,  30.21001000, 17.90520000);
  mol.emplace_back(AtomicNumber(1 ), 17.78336000,  30.05335000, 18.89396000);
  mol.emplace_back(AtomicNumber(1 ), 16.72586000,  28.03197000, 18.38075000);
  mol.emplace_back(AtomicNumber(1 ), 18.21258000,  28.00240000, 17.26090000);
  mol.emplace_back(AtomicNumber(1 ), 18.21242000,  27.00786000, 18.83485000);
  mol.emplace_back(AtomicNumber(1 ), 19.66208000,  30.56192000, 20.64032000);
  mol.emplace_back(AtomicNumber(1 ), 20.03907000,  30.81855000, 22.84239000);
  mol.emplace_back(AtomicNumber(1 ), 18.35667000,  29.92134000, 23.11588000);
  mol.emplace_back(AtomicNumber(1 ), 20.84605000,  28.11284000, 21.51183000);
  mol.emplace_back(AtomicNumber(1 ), 21.27080000,  26.67093000, 24.23118000);
  mol.emplace_back(AtomicNumber(1 ), 21.52124000,  25.66819000, 21.29584000);
  mol.emplace_back(AtomicNumber(1 ), 21.84640000,  24.69685000, 22.89930000);
  mol.emplace_back(AtomicNumber(1 ), 19.48390000,  24.78208000, 23.52050000);
  mol.emplace_back(AtomicNumber(1 ), 19.06355000,  25.92146000, 22.14091000);
  mol.emplace_back(AtomicNumber(1 ), 20.67366000,  23.56776000, 21.18740000);
  mol.emplace_back(AtomicNumber(1 ), 19.12328000,  23.04012000, 22.10104000);
  mol.emplace_back(AtomicNumber(1 ), 19.27583000,  25.06568000, 19.72050000);
  mol.emplace_back(AtomicNumber(1 ), 18.69614000,  23.32147000, 19.61913000);
  mol.emplace_back(AtomicNumber(1 ), 16.82157000,  23.75135000, 21.06229000);
  mol.emplace_back(AtomicNumber(1 ), 17.57303000,  25.26808000, 21.75481000);
  mol.emplace_back(AtomicNumber(1 ), 23.50379000,  26.32231000, 24.78186000);
  mol.emplace_back(AtomicNumber(1 ), 25.62570000,  27.80195000, 23.04448000);
  mol.emplace_back(AtomicNumber(1 ), 25.59280000,  27.86920000, 25.55934000);
  mol.emplace_back(AtomicNumber(1 ), 25.99581000,  26.05279000, 25.55068000);
  mol.emplace_back(AtomicNumber(1 ), 28.04340000,  26.48947000, 24.13981000);
  mol.emplace_back(AtomicNumber(1 ), 27.62638000,  28.30048000, 24.06083000);
  mol.emplace_back(AtomicNumber(1 ), 26.94550000,  29.16713000, 26.36278000);
  mol.emplace_back(AtomicNumber(1 ), 28.21391000,  28.82379000, 27.70020000);
  mol.emplace_back(AtomicNumber(1 ), 26.80853000,  26.66851000, 21.52583000);
  mol.emplace_back(AtomicNumber(1 ), 25.98879000,  23.74829000, 20.97202000);
  mol.emplace_back(AtomicNumber(1 ), 28.00304000,  25.60632000, 19.54789000);
  mol.emplace_back(AtomicNumber(1 ), 27.07551000,  24.04082000, 18.94346000);
  mol.emplace_back(AtomicNumber(1 ), 25.92228000,  26.82622000, 19.53759000);
  mol.emplace_back(AtomicNumber(1 ), 26.30075000,  25.08094000, 16.93627000);
  mol.emplace_back(AtomicNumber(1 ), 25.33975000,  26.67319000, 17.01705000);
  mol.emplace_back(AtomicNumber(1 ), 27.17121000,  26.67291000, 17.35197000);
  mol.emplace_back(AtomicNumber(1 ), 24.69951000,  24.05665000, 19.06814000);
  mol.emplace_back(AtomicNumber(1 ), 24.08608000,  25.50966000, 20.05696000);
  mol.emplace_back(AtomicNumber(1 ), 23.92199000,  25.50951000, 18.20238000);
  mol.emplace_back(AtomicNumber(1 ), 27.10774000,  21.94139000, 21.17985000);
  mol.emplace_back(AtomicNumber(1 ), 29.63959000,  22.30342000, 22.98916000);
  mol.emplace_back(AtomicNumber(1 ), 27.75175000,  21.02419000, 23.90841000);
  mol.emplace_back(AtomicNumber(1 ), 27.73852000,  19.97658000, 22.28454000);
  mol.emplace_back(AtomicNumber(1 ), 29.26962000,  18.54186000, 23.02423000);
  mol.emplace_back(AtomicNumber(1 ), 30.29554000,  19.88429000, 23.83435000);
  mol.emplace_back(AtomicNumber(1 ), 27.49987000,  20.61753000, 24.69411000);
  mol.emplace_back(AtomicNumber(1 ), 31.56853000,  21.41977000, 22.43392000);
  mol.emplace_back(AtomicNumber(1 ), 32.13736000,  21.42139000, 19.66578000);
  mol.emplace_back(AtomicNumber(1 ), 33.51443000,  20.28507000, 22.29016000);
  mol.emplace_back(AtomicNumber(1 ), 34.38081000,  20.10510000, 20.62992000);
  mol.emplace_back(AtomicNumber(1 ), 35.77759000,  21.50109000, 22.38332000);
  mol.emplace_back(AtomicNumber(1 ), 31.41754000,  18.44239000, 21.44352000);
  mol.emplace_back(AtomicNumber(1 ), 32.45311000,  16.79775000, 19.12557000);
  mol.emplace_back(AtomicNumber(1 ), 31.53495000,  16.24771000, 20.71286000);
  mol.emplace_back(AtomicNumber(1 ), 29.44788000,  18.46858000, 19.70471000);
  mol.emplace_back(AtomicNumber(1 ), 27.83584000,  16.33254000, 18.24038000);
  mol.emplace_back(AtomicNumber(1 ), 27.15703000,  19.30694000, 18.86194000);
  mol.emplace_back(AtomicNumber(1 ), 26.00557000,  17.91468000, 18.27877000);
  mol.emplace_back(AtomicNumber(1 ), 27.57831000,  17.79221000, 20.97095000);
  mol.emplace_back(AtomicNumber(1 ), 25.90957000,  18.53852000, 20.84856000);
  mol.emplace_back(AtomicNumber(1 ), 25.00681000,  16.33466000, 20.82034000);
  mol.emplace_back(AtomicNumber(1 ), 25.81228000,  16.04053000, 19.20671000);
  mol.emplace_back(AtomicNumber(1 ), 27.23964000,  15.89540000, 21.92261000);
  mol.emplace_back(AtomicNumber(1 ), 26.80521000,  14.07495000, 18.76183000);
  mol.emplace_back(AtomicNumber(1 ), 28.05107000,  12.82504000, 19.39506000);
  mol.emplace_back(AtomicNumber(1 ), 28.70250000,  14.51665000, 22.46434000);
  mol.emplace_back(AtomicNumber(1 ), 26.69749000,  17.02347000, 16.21612000);
  mol.emplace_back(AtomicNumber(1 ), 28.61455000,  18.60344000, 14.34488000);
  mol.emplace_back(AtomicNumber(1 ), 27.44147000,  17.45637000, 12.47555000);
  mol.emplace_back(AtomicNumber(1 ), 25.88429000,  15.51181000, 13.92426000);
  mol.emplace_back(AtomicNumber(1 ), 27.62662000,  14.90075000, 13.57822000);
  mol.emplace_back(AtomicNumber(1 ), 28.89331000,  15.87818000, 12.62669000);
  mol.emplace_back(AtomicNumber(1 ), 28.89300000,  15.91724000, 14.48811000);
  mol.emplace_back(AtomicNumber(1 ), 27.42868000,  19.93314000, 12.71231000);
  mol.emplace_back(AtomicNumber(1 ), 25.47383000,  21.85707000, 13.92380000);
  mol.emplace_back(AtomicNumber(1 ), 26.15186000,  21.14181000, 10.97382000);
  mol.emplace_back(AtomicNumber(1 ), 25.32161000,  22.72760000, 11.71712000);
  mol.emplace_back(AtomicNumber(1 ), 28.21409000,  21.82735000, 12.13980000);
  mol.emplace_back(AtomicNumber(1 ), 26.73522000,  23.61847000, 10.07129000);
  mol.emplace_back(AtomicNumber(1 ), 28.31786000,  24.28752000, 10.78747000);
  mol.emplace_back(AtomicNumber(1 ), 28.31757000,  22.65442000,  9.89337000);
  mol.emplace_back(AtomicNumber(1 ), 26.29441000,  23.97293000, 13.19507000);
  mol.emplace_back(AtomicNumber(1 ), 27.90581000,  23.36546000, 13.90195000);
  mol.emplace_back(AtomicNumber(1 ), 27.90556000,  24.64974000, 12.55397000);
  mol.emplace_back(AtomicNumber(1 ), 25.16148000,  18.76108000, 11.92983000);
  mol.emplace_back(AtomicNumber(1 ), 22.27785000,  19.34196000, 11.47176000);
  mol.emplace_back(AtomicNumber(1 ), 22.03031000,  16.85439000, 11.01622000);
  mol.emplace_back(AtomicNumber(1 ), 23.49620000,  17.62571000, 10.14473000);
  mol.emplace_back(AtomicNumber(1 ), 23.31736000,  15.61859000, 12.25620000);
  mol.emplace_back(AtomicNumber(1 ), 24.07292000,  18.00815000, 14.25067000);
  mol.emplace_back(AtomicNumber(1 ), 21.82746000,  16.67313000, 15.55200000);
  mol.emplace_back(AtomicNumber(1 ), 24.24350000,  18.15906000, 16.70079000);
  mol.emplace_back(AtomicNumber(1 ), 22.88468000,  17.26287000, 17.72223000);
  mol.emplace_back(AtomicNumber(1 ), 25.73151000,  16.71054000, 17.68179000);
  mol.emplace_back(AtomicNumber(1 ), 22.79520000,  20.08887000, 15.40237000);
  mol.emplace_back(AtomicNumber(1 ), 20.40934000,  20.91564000, 17.14358000);
  mol.emplace_back(AtomicNumber(1 ), 22.76684000,  22.40123000, 15.67932000);
  mol.emplace_back(AtomicNumber(1 ), 21.48316000,  23.19014000, 16.82802000);
  mol.emplace_back(AtomicNumber(1 ), 24.70663000,  21.42630000, 16.63874000);
  mol.emplace_back(AtomicNumber(1 ), 21.40589000,  22.52568000, 19.18414000);
  mol.emplace_back(AtomicNumber(1 ), 26.09018000,  20.93441000, 18.66387000);
  mol.emplace_back(AtomicNumber(1 ), 22.78241000,  22.05144000, 21.19932000);
  mol.emplace_back(AtomicNumber(1 ), 26.27135000,  21.01733000, 21.13854000);
  mol.emplace_back(AtomicNumber(1 ), 20.54780000,  19.96437000, 13.91382000);
  mol.emplace_back(AtomicNumber(1 ), 19.16541000,  20.64089000, 12.11365000);
  mol.emplace_back(AtomicNumber(1 ), 17.44844000,  20.33436000, 14.07958000);
  mol.emplace_back(AtomicNumber(1 ), 17.38271000,  22.23097000, 14.05910000);
  mol.emplace_back(AtomicNumber(1 ), 17.88482000,  20.08983000, 11.16648000);
  mol.emplace_back(AtomicNumber(1 ), 16.24550000,  20.53701000, 10.37395000);
  mol.emplace_back(AtomicNumber(1 ), 21.44030000,  22.18701000, 12.47195000);
  mol.emplace_back(AtomicNumber(1 ), 20.55087000,  25.15026000, 12.02394000);
  mol.emplace_back(AtomicNumber(1 ), 23.22613000,  23.55807000, 11.94828000);
  mol.emplace_back(AtomicNumber(1 ), 22.08664000,  25.78295000, 13.89738000);
  mol.emplace_back(AtomicNumber(1 ), 22.17064000,  23.95450000, 14.25649000);
  mol.emplace_back(AtomicNumber(1 ), 22.60653000,  26.34191000, 10.93801000);
  mol.emplace_back(AtomicNumber(1 ), 24.02899000,  26.23456000, 12.13400000);
  mol.emplace_back(AtomicNumber(1 ), 24.02868000,  25.19201000, 10.59143000);
  mol.emplace_back(AtomicNumber(1 ), 24.58763000,  23.95200000, 14.35274000);
  mol.emplace_back(AtomicNumber(1 ), 24.72682000,  25.55522000, 13.41701000);
  mol.emplace_back(AtomicNumber(1 ), 24.12814000,  25.55495000, 15.17996000);
  mol.emplace_back(AtomicNumber(1 ), 19.73242000,  25.62464000, 10.05274000);
  mol.emplace_back(AtomicNumber(1 ), 20.51742000,  23.80024000,  7.62943000);
  mol.emplace_back(AtomicNumber(1 ), 18.10345000,  25.50951000,  8.40206000);
  mol.emplace_back(AtomicNumber(1 ), 18.55720000,  24.92100000,  6.60289000);
  mol.emplace_back(AtomicNumber(1 ), 18.43832000,  22.56229000,  7.50483000);
  mol.emplace_back(AtomicNumber(1 ), 17.59468000,  23.26386000,  9.00469000);
  mol.emplace_back(AtomicNumber(1 ), 17.45429000,  24.79490000,  5.72239000);
  mol.emplace_back(AtomicNumber(1 ), 15.64785000,  24.40582000,  5.40491000);
  mol.emplace_back(AtomicNumber(1 ), 20.43959000,  25.04973000,  5.47935000);
  mol.emplace_back(AtomicNumber(1 ), 22.74776000,  26.72543000,  5.54477000);
  mol.emplace_back(AtomicNumber(1 ), 22.07794000,  27.56435000,  3.23189000);
  mol.emplace_back(AtomicNumber(1 ), 22.23361000,  25.70098000,  3.50529000);
  mol.emplace_back(AtomicNumber(1 ), 19.84947000,  25.37356000,  3.43631000);
  mol.emplace_back(AtomicNumber(1 ), 19.48567000,  27.17711000,  3.62916000);
  mol.emplace_back(AtomicNumber(1 ), 20.17235000,  27.69994000,  1.41879000);
  mol.emplace_back(AtomicNumber(1 ), 21.26104000,  26.13658000,  1.32616000);
  mol.emplace_back(AtomicNumber(1 ), 19.47140000,  25.84523000, -0.21665000);
  mol.emplace_back(AtomicNumber(1 ), 19.12683000,  24.78774000,  1.27912000);
  mol.emplace_back(AtomicNumber(1 ), 17.07534000,  25.77737000,  1.18789000);
  mol.emplace_back(AtomicNumber(1 ), 17.92191000,  27.18442000,  1.99305000);
  mol.emplace_back(AtomicNumber(1 ), 23.10931000,  28.82877000,  5.40658000);
  mol.emplace_back(AtomicNumber(1 ), 22.96667000,  30.98203000,  5.84756000);
  mol.emplace_back(AtomicNumber(1 ), 20.10224000,  30.57837000,  4.59359000);
  mol.emplace_back(AtomicNumber(1 ), 20.68585000,  32.25700000,  5.12571000);
  mol.emplace_back(AtomicNumber(1 ), 22.04590000,  32.44485000,  3.25234000);
  mol.emplace_back(AtomicNumber(1 ), 22.84443000,  30.78661000,  3.62892000);
  mol.emplace_back(AtomicNumber(1 ), 22.99551000,  30.14889000,  1.75015000);
  mol.emplace_back(AtomicNumber(1 ), 22.05273000,  29.06852000,  7.88079000);
  mol.emplace_back(AtomicNumber(1 ), 20.36482000,  30.65306000,  9.78809000);
  mol.emplace_back(AtomicNumber(1 ), 22.61260000,  28.44849000, 10.00279000);
  mol.emplace_back(AtomicNumber(1 ), 21.69371000,  29.10612000, 11.47677000);
  mol.emplace_back(AtomicNumber(1 ), 20.87645000,  27.00226000, 10.14842000);
  mol.emplace_back(AtomicNumber(1 ), 21.26702000,  31.74688000, 11.69626000);
  mol.emplace_back(AtomicNumber(1 ), 24.15722000,  32.77723000, 11.25387000);
  mol.emplace_back(AtomicNumber(1 ), 21.73218000,  34.37293000, 12.57967000);
  mol.emplace_back(AtomicNumber(1 ), 21.41723000,  34.75381000, 10.37842000);
  mol.emplace_back(AtomicNumber(1 ), 24.47499000,  34.45828000, 13.25825000);
  mol.emplace_back(AtomicNumber(1 ), 24.27788000,  35.78910000, 11.97164000);
  mol.emplace_back(AtomicNumber(1 ), 23.19306000,  35.78878000, 13.48478000);
  mol.emplace_back(AtomicNumber(1 ), 25.56512000,  32.16126000, 12.87549000);
  mol.emplace_back(AtomicNumber(1 ), 24.24360000,  31.14411000, 15.49354000);
  mol.emplace_back(AtomicNumber(1 ), 27.26902000,  31.18517000, 14.79264000);
  mol.emplace_back(AtomicNumber(1 ), 26.21763000,  30.05904000, 15.96944000);
  mol.emplace_back(AtomicNumber(1 ), 26.49007000,  29.94221000, 12.93729000);
  mol.emplace_back(AtomicNumber(1 ), 26.96139000,  27.39213000, 14.12374000);
  mol.emplace_back(AtomicNumber(1 ), 28.23827000,  28.58237000, 13.47706000);
  mol.emplace_back(AtomicNumber(1 ), 27.79684000,  28.58230000, 15.28580000);
  mol.emplace_back(AtomicNumber(1 ), 25.08153000,  27.67557000, 13.95283000);
  mol.emplace_back(AtomicNumber(1 ), 24.35654000,  29.11822000, 14.87935000);
  mol.emplace_back(AtomicNumber(1 ), 24.36786000,  29.11807000, 13.01755000);
  mol.emplace_back(AtomicNumber(1 ), 25.15304000,  32.14746000, 17.51280000);
  mol.emplace_back(AtomicNumber(1 ), 26.33772000,  34.99577000, 16.92913000);
  mol.emplace_back(AtomicNumber(1 ), 24.94157000,  33.66831000, 19.42909000);
  mol.emplace_back(AtomicNumber(1 ), 25.70423000,  35.34932000, 19.44642000);
  mol.emplace_back(AtomicNumber(1 ), 24.59888000,  37.19735000, 17.95337000);
  mol.emplace_back(AtomicNumber(1 ), 22.56016000,  33.45929000, 17.95884000);
  mol.emplace_back(AtomicNumber(1 ), 22.28838000,  37.58310000, 16.81209000);
  mol.emplace_back(AtomicNumber(1 ), 28.43482000,  35.35154000, 17.60122000);
  mol.emplace_back(AtomicNumber(1 ), 29.78747000,  33.07876000, 19.23218000);
  mol.emplace_back(AtomicNumber(1 ), 30.79328000,  33.49203000, 16.96373000);
  mol.emplace_back(AtomicNumber(1 ), 30.82802000,  35.34631000, 17.30158000);
  mol.emplace_back(AtomicNumber(1 ), 32.47608000,  34.91217000, 19.22445000);
  mol.emplace_back(AtomicNumber(1 ), 31.54398000,  32.28832000, 19.05091000);
  mol.emplace_back(AtomicNumber(1 ), 33.14023000,  32.18045000, 18.09926000);
  mol.emplace_back(AtomicNumber(1 ), 33.13995000,  32.82006000, 19.84778000);
  mol.emplace_back(AtomicNumber(1 ), 33.75743000,  33.51435000, 16.81072000);
  mol.emplace_back(AtomicNumber(1 ), 32.85190000,  35.09368000, 16.42208000);
  mol.emplace_back(AtomicNumber(1 ), 34.22688000,  35.09339000, 17.67740000);
  mol.emplace_back(AtomicNumber(1 ), 29.92241000,  33.47837000, 21.32857000);
  mol.emplace_back(AtomicNumber(1 ), 30.44952000,  36.49070000, 22.17268000);
  mol.emplace_back(AtomicNumber(1 ), 29.54177000,  34.07349000, 23.75907000);
  mol.emplace_back(AtomicNumber(1 ), 28.86188000,  36.27883000, 25.29047000);
  mol.emplace_back(AtomicNumber(1 ), 30.43112000,  35.29219000, 25.46169000);
  mol.emplace_back(AtomicNumber(1 ), 30.43083000,  36.87983000, 24.48917000);
  mol.emplace_back(AtomicNumber(1 ), 27.50450000,  34.46586000, 22.90761000);
  mol.emplace_back(AtomicNumber(1 ), 27.51667000,  36.05130000, 23.88307000);
  mol.emplace_back(AtomicNumber(1 ), 28.05191000,  36.05101000, 22.09983000);
  mol.emplace_back(AtomicNumber(1 ), 32.45863000,  36.86975000, 23.10290000);
  mol.emplace_back(AtomicNumber(1 ), 34.44353000,  34.45926000, 23.05107000);
  mol.emplace_back(AtomicNumber(1 ), 34.50469000,  37.45367000, 22.53709000);
  mol.emplace_back(AtomicNumber(1 ), 35.82839000,  36.90206000, 23.72856000);
  mol.emplace_back(AtomicNumber(1 ), 36.79629000,  35.34971000, 22.11935000);
  mol.emplace_back(AtomicNumber(1 ), 35.24520000,  33.99085000, 21.13438000);
  mol.emplace_back(AtomicNumber(1 ), 33.91942000,  35.29575000, 21.20393000);
  mol.emplace_back(AtomicNumber(1 ), 35.15367000,  35.29565000, 19.81000000);
  mol.emplace_back(AtomicNumber(1 ), 35.70288000,  37.71506000, 20.45219000);
  mol.emplace_back(AtomicNumber(1 ), 37.17800000,  37.74604000, 21.58725000);
  mol.emplace_back(AtomicNumber(1 ), 37.17770000,  36.63843000, 20.09072000);
  mol.emplace_back(AtomicNumber(1 ), 35.36861000,  33.68078000, 24.75464000);
  mol.emplace_back(AtomicNumber(1 ), 34.58186000,  34.92166000, 27.53255000);
  mol.emplace_back(AtomicNumber(1 ), 33.73517000,  32.82842000, 27.86126000);
  mol.emplace_back(AtomicNumber(1 ), 34.59356000,  32.15384000, 26.32843000);
  mol.emplace_back(AtomicNumber(1 ), 36.79769000,  32.35923000, 27.82388000);
  mol.emplace_back(AtomicNumber(1 ), 35.58453000,  32.35593000, 29.24646000);
  mol.emplace_back(AtomicNumber(1 ), 36.63610000,  30.08684000, 28.14678000);
  mol.emplace_back(AtomicNumber(1 ), 34.92697000,  30.18928000, 28.86863000);
  mol.emplace_back(AtomicNumber(1 ), 35.17467000,  30.89941000, 25.87296000);
  mol.emplace_back(AtomicNumber(1 ), 34.48311000,  28.81738000, 28.71029000);
  mol.emplace_back(AtomicNumber(1 ), 33.55008000,  27.51093000, 27.74175000);
  mol.emplace_back(AtomicNumber(1 ), 34.31025000,  29.20874000, 24.57404000);
  mol.emplace_back(AtomicNumber(1 ), 36.56438000,  36.24375000, 27.74670000);
  mol.emplace_back(AtomicNumber(1 ), 39.24161000,  34.86820000, 26.92576000);
  mol.emplace_back(AtomicNumber(1 ), 38.19290000,  37.59985000, 27.68630000);
  mol.emplace_back(AtomicNumber(1 ), 40.03139000,  37.07068000, 28.02053000);
  mol.emplace_back(AtomicNumber(1 ), 40.58679000,  37.02079000, 25.85733000);
  mol.emplace_back(AtomicNumber(1 ), 38.77513000,  35.51702000, 25.04321000);
  mol.emplace_back(AtomicNumber(1 ), 37.59002000,  36.95240000, 25.02071000);
  mol.emplace_back(AtomicNumber(1 ), 39.07386000,  36.95226000, 23.89615000);
  mol.emplace_back(AtomicNumber(1 ), 38.22761000,  39.11527000, 25.83784000);
  mol.emplace_back(AtomicNumber(1 ), 39.83038000,  39.37349000, 26.74874000);
  mol.emplace_back(AtomicNumber(1 ), 39.83014000,  39.20185000, 24.89484000);
  mol.emplace_back(AtomicNumber(1 ), 40.88378000,  34.43713000, 28.13762000);
  mol.emplace_back(AtomicNumber(1 ), 40.02948000,  33.49563000, 30.95600000);
  mol.emplace_back(AtomicNumber(1 ), 41.22197000,  32.09393000, 28.96058000);
  mol.emplace_back(AtomicNumber(1 ), 42.74383000,  32.84181000, 29.83653000);
  mol.emplace_back(AtomicNumber(1 ), 42.38265000,  31.69109000, 31.78804000);
  mol.emplace_back(AtomicNumber(1 ), 40.54825000,  31.42532000, 31.40874000);
  mol.emplace_back(AtomicNumber(1 ), 41.28651000,  29.31861000, 31.04627000);
  mol.emplace_back(AtomicNumber(1 ), 41.61483000,  30.06860000, 29.35128000);
  mol.emplace_back(AtomicNumber(1 ), 43.67668000,  29.32529000, 31.49586000);
  mol.emplace_back(AtomicNumber(1 ), 42.70281000,  30.50960000, 28.18936000);
  mol.emplace_back(AtomicNumber(1 ), 44.46415000,  30.44488000, 27.54987000);
  mol.emplace_back(AtomicNumber(1 ), 45.83392000,  29.44432000, 30.74538000);
  mol.emplace_back(AtomicNumber(1 ), 40.11945000,  35.33656000, 32.10003000);
  mol.emplace_back(AtomicNumber(1 ), 41.81486000,  36.22932000, 33.77955000);
  mol.emplace_back(AtomicNumber(1 ), 42.90134000,  36.62534000, 32.26401000);
  mol.emplace_back(AtomicNumber(1 ), 41.13110000,  37.87717000, 34.64228000);
  mol.emplace_back(AtomicNumber(1 ), 41.13053000,  40.60310000, 33.62542000);
  mol.emplace_back(AtomicNumber(1 ), 39.43295000,  39.92894000, 33.30960000);
  mol.emplace_back(AtomicNumber(1 ), 41.76077000,  39.11904000, 35.91137000);

  // Convert to Bohr
  for( auto& atom : mol ) {
    atom.x *= 1.8897161646321;
    atom.y *= 1.8897161646321;
    atom.z *= 1.8897161646321;
  }
#else
  mol.emplace_back(AtomicNumber(7), 5.1665108539259997e+01, 4.6166005911269998e+01, 4.9397437352460001e+00);
  mol.emplace_back(AtomicNumber(6), 4.9635542827073998e+01, 4.8023606558456997e+01, 5.3706012607379998e+00);
  mol.emplace_back(AtomicNumber(6), 5.0858195541957002e+01, 5.0340410620970999e+01, 6.6726224671590000e+00);
  mol.emplace_back(AtomicNumber(8), 5.2696898929253997e+01, 5.0007818846907000e+01, 8.0559018911070002e+00);
  mol.emplace_back(AtomicNumber(6), 4.7454799035767998e+01, 4.7016382606320001e+01, 6.8956101338610001e+00);
  mol.emplace_back(AtomicNumber(6), 4.7910222999117003e+01, 4.6978588086540000e+01, 9.7018532275260014e+00);
  mol.emplace_back(AtomicNumber(16), 4.5221142916769999e+01, 4.5275944970451000e+01, 1.1156942239056001e+01);
  mol.emplace_back(AtomicNumber(6), 4.6198131253082998e+01, 4.5323188120176006e+01, 1.4399712036180000e+01);
  mol.emplace_back(AtomicNumber(7), 4.9765933920315000e+01, 5.2477690714529999e+01, 6.1567272721619997e+00);
  mol.emplace_back(AtomicNumber(6), 5.0739142804650001e+01, 5.4841737926769000e+01, 7.3661519051219999e+00);
  mol.emplace_back(AtomicNumber(6), 4.9321848312900002e+01, 5.5280154356216997e+01, 9.8303545947779991e+00);
  mol.emplace_back(AtomicNumber(8), 4.6988036716484999e+01, 5.4847407104736000e+01, 1.0072239521369999e+01);
  mol.emplace_back(AtomicNumber(6), 5.0518044863937000e+01, 5.6971459116371996e+01, 5.4896539980449992e+00);
  mol.emplace_back(AtomicNumber(6), 5.0799614036298003e+01, 5.9613296048994002e+01, 6.4420758965009997e+00);
  mol.emplace_back(AtomicNumber(6), 5.0618200341354004e+01, 6.1533257653817998e+01, 4.2896779950299999e+00);
  mol.emplace_back(AtomicNumber(8), 5.2502257152387003e+01, 6.2663313795239993e+01, 3.5337875994300001e+00);
  mol.emplace_back(AtomicNumber(7), 4.8305175730818000e+01, 6.1856400797936992e+01, 3.4128451361340000e+00);
  mol.emplace_back(AtomicNumber(7), 5.0737253078660999e+01, 5.6041713929783995e+01, 1.1748426473613000e+01);
  mol.emplace_back(AtomicNumber(6), 4.9576961321414998e+01, 5.6801383777361998e+01, 1.4167275739533000e+01);
  mol.emplace_back(AtomicNumber(6), 5.0799614036298003e+01, 5.9390308382291998e+01, 1.4857025725518000e+01);
  mol.emplace_back(AtomicNumber(8), 5.2734693449033998e+01, 5.9925100837178995e+01, 1.3726969584096000e+01);
  mol.emplace_back(AtomicNumber(6), 4.9782941454216001e+01, 5.4896539980450001e+01, 1.6336681174904999e+01);
  mol.emplace_back(AtomicNumber(6), 5.2553279754089999e+01, 5.4325842731771999e+01, 1.7005644175011000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8171005185599000e+01, 5.2479580440519001e+01, 1.5660159270843002e+01);
  mol.emplace_back(AtomicNumber(6), 5.2849966734363001e+01, 5.3076733853043002e+01, 1.9685275627412999e+01);
  mol.emplace_back(AtomicNumber(7), 4.9537277075645996e+01, 6.0654535068933001e+01, 1.6574786649519002e+01);
  mol.emplace_back(AtomicNumber(6), 5.0591744177507998e+01, 6.3184878168204001e+01, 1.7379809920832997e+01);
  mol.emplace_back(AtomicNumber(6), 5.1307950327339000e+01, 6.3045038445018001e+01, 2.0125581782850002e+01);
  mol.emplace_back(AtomicNumber(8), 4.9794279810150002e+01, 6.1941438467441998e+01, 2.1533427644654999e+01);
  mol.emplace_back(AtomicNumber(6), 4.8556509287354999e+01, 6.5191767168521991e+01, 1.6905488697593999e+01);
  mol.emplace_back(AtomicNumber(6), 4.7787390809831997e+01, 6.5401526753300999e+01, 1.4171055191511000e+01);
  mol.emplace_back(AtomicNumber(6), 4.5631213456382994e+01, 6.4186432942373997e+01, 1.3299891510582000e+01);
  mol.emplace_back(AtomicNumber(6), 4.9389878448504000e+01, 6.6794254807193994e+01, 1.2547780566960000e+01);
  mol.emplace_back(AtomicNumber(6), 4.4998155250068002e+01, 6.4309265131659004e+01, 1.0727974439553000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8773827776089995e+01, 6.6881182202688009e+01, 9.9531867840629999e+00);
  mol.emplace_back(AtomicNumber(6), 4.6525053849180004e+01, 6.5720890445441995e+01, 9.1708402246169989e+00);
  mol.emplace_back(AtomicNumber(7), 5.3403656449140001e+01, 6.4142969244626997e+01, 2.0968399573944001e+01);
  mol.emplace_back(AtomicNumber(6), 5.4055611915345004e+01, 6.4184543216385009e+01, 2.3627244040467001e+01);
  mol.emplace_back(AtomicNumber(6), 5.4117972872982001e+01, 6.7011573295928997e+01, 2.4377465258099999e+01);
  mol.emplace_back(AtomicNumber(8), 5.5788490647257994e+01, 6.8224777380866996e+01, 2.3281424184479999e+01);
  mol.emplace_back(AtomicNumber(6), 5.6621859808407002e+01, 6.2960000775513002e+01, 2.4214948823046001e+01);
  mol.emplace_back(AtomicNumber(6), 5.7090511853678997e+01, 6.3105509676665996e+01, 2.7030640546655999e+01);
  mol.emplace_back(AtomicNumber(6), 5.6610521452473002e+01, 6.0165096037782000e+01, 2.3341895416128001e+01);
  mol.emplace_back(AtomicNumber(7), 5.2441785920739001e+01, 6.7778802047462989e+01, 2.5964835088859999e+01);
  mol.emplace_back(AtomicNumber(6), 5.2328402361399000e+01, 7.0515125279534999e+01, 2.6726394662427001e+01);
  mol.emplace_back(AtomicNumber(6), 5.3798609180841005e+01, 7.0817481437775001e+01, 2.9139574750380000e+01);
  mol.emplace_back(AtomicNumber(8), 5.3314839327657005e+01, 6.9453099273717001e+01, 3.1012293205479004e+01);
  mol.emplace_back(AtomicNumber(6), 4.9546725705591001e+01, 7.1212434169475998e+01, 2.7036309724622999e+01);
  mol.emplace_back(AtomicNumber(6), 4.8913667499276002e+01, 7.3961985483471011e+01, 2.7618345329235002e+01);
  mol.emplace_back(AtomicNumber(6), 4.6011048380171999e+01, 7.4258672463744006e+01, 2.7669367930937998e+01);
  mol.emplace_back(AtomicNumber(6), 4.5098310727485000e+01, 7.6955311450046992e+01, 2.7871568611761003e+01);
  mol.emplace_back(AtomicNumber(7), 4.2282619003874998e+01, 7.6949642272079998e+01, 2.8170145318023000e+01);
  mol.emplace_back(AtomicNumber(7), 5.5607076952313996e+01, 7.2622169757270001e+01, 2.9188707626094001e+01);
  mol.emplace_back(AtomicNumber(6), 5.7116968017525004e+01, 7.3024681392927008e+01, 3.1486614428717999e+01);
  mol.emplace_back(AtomicNumber(6), 5.6056831737696001e+01, 7.5284793675770999e+01, 3.2945482892226003e+01);
  mol.emplace_back(AtomicNumber(8), 5.4518594782650005e+01, 7.6656734743784995e+01, 3.1858890448551005e+01);
  mol.emplace_back(AtomicNumber(6), 5.9987461794815999e+01, 7.3470656726331001e+01, 3.0800643894710998e+01);
  mol.emplace_back(AtomicNumber(8), 5.9974233712892996e+01, 7.6074699139172992e+01, 2.9903024049936000e+01);
  mol.emplace_back(AtomicNumber(6), 6.0962560405139996e+01, 7.1751006076341000e+01, 2.8669032979118999e+01);
  mol.emplace_back(AtomicNumber(7), 5.6941223500548006e+01, 7.5719430653241005e+01, 3.5228271886937996e+01);
  mol.emplace_back(AtomicNumber(6), 5.5949117356323001e+01, 7.7818916227019997e+01, 3.6787295827862998e+01);
  mol.emplace_back(AtomicNumber(6), 5.6833509119174998e+01, 8.0385164120081996e+01, 3.5874558175176006e+01);
  mol.emplace_back(AtomicNumber(8), 5.5909433110553998e+01, 8.2335361340730003e+01, 3.6817531443687002e+01);
  mol.emplace_back(AtomicNumber(6), 5.6538711864890999e+01, 7.7270895690209997e+01, 3.9567082757681995e+01);
  mol.emplace_back(AtomicNumber(6), 5.5147873536986999e+01, 7.5063695735058005e+01, 4.0782176568608996e+01);
  mol.emplace_back(AtomicNumber(6), 5.5384089285611999e+01, 7.5116608062750004e+01, 4.3643221715955001e+01);
  mol.emplace_back(AtomicNumber(6), 5.2345409895300001e+01, 7.5061806009068988e+01, 4.0115103294492002e+01);
  mol.emplace_back(AtomicNumber(7), 5.8564498125099000e+01, 8.0447525077718993e+01, 3.4011288350021999e+01);
  mol.emplace_back(AtomicNumber(6), 5.9378970026358004e+01, 8.3034559956659990e+01, 3.3170360284917003e+01);
  mol.emplace_back(AtomicNumber(6), 5.8118522791695000e+01, 8.3811237338138994e+01, 3.0759069922953003e+01);
  mol.emplace_back(AtomicNumber(8), 5.8972678938723000e+01, 8.5544116070051999e+01, 2.9415474744774002e+01);
  mol.emplace_back(AtomicNumber(6), 6.2321273391230996e+01, 8.2992985984901992e+01, 3.2966269878105003e+01);
  mol.emplace_back(AtomicNumber(8), 6.2689769959086000e+01, 8.1384829168262996e+01, 3.0736393211085002e+01);
  mol.emplace_back(AtomicNumber(6), 6.3602507611772992e+01, 8.1861040117491001e+01, 3.5284963666608000e+01);
  mol.emplace_back(AtomicNumber(7), 5.6164546119069001e+01, 8.2530003117597005e+01, 3.0018297335265000e+01);
  mol.emplace_back(AtomicNumber(6), 5.4760479709242006e+01, 8.3072354476439997e+01, 2.7737398066542003e+01);
  mol.emplace_back(AtomicNumber(6), 5.5943448178356000e+01, 8.2216308603423002e+01, 2.5309100170677002e+01);
  mol.emplace_back(AtomicNumber(8), 5.5215903672591004e+01, 8.3112038722209007e+01, 2.3245519390689001e+01);
  mol.emplace_back(AtomicNumber(7), 5.7755695401806996e+01, 8.0545790829146995e+01, 2.5501852221554998e+01);
  mol.emplace_back(AtomicNumber(6), 5.8942443322898995e+01, 7.9391168249868002e+01, 2.3302211170358998e+01);
  mol.emplace_back(AtomicNumber(6), 5.7559163898950999e+01, 7.6847597068673991e+01, 2.2922376246570000e+01);
  mol.emplace_back(AtomicNumber(8), 5.7169880345217003e+01, 7.5572032026098995e+01, 2.4817771413536999e+01);
  mol.emplace_back(AtomicNumber(6), 6.1741127512607996e+01, 7.8833699083113004e+01, 2.3631023492445003e+01);
  mol.emplace_back(AtomicNumber(6), 6.2890080913920002e+01, 7.7641281984054004e+01, 2.1215953678503002e+01);
  mol.emplace_back(AtomicNumber(6), 6.5690654829617998e+01, 7.7098930625210997e+01, 2.1675157093830002e+01);
  mol.emplace_back(AtomicNumber(6), 6.7300701372245996e+01, 7.7189637472683003e+01, 1.9350794127360000e+01);
  mol.emplace_back(AtomicNumber(7), 6.6329382213900004e+01, 7.5726989557197001e+01, 1.7198396225889002e+01);
  mol.emplace_back(AtomicNumber(7), 5.6999805006206998e+01, 7.6227766944281996e+01, 2.0571557116253999e+01);
  mol.emplace_back(AtomicNumber(6), 5.5826285167038002e+01, 7.3737108090780012e+01, 2.0131250960817002e+01);
  mol.emplace_back(AtomicNumber(6), 5.7625304308566001e+01, 7.2302806065129005e+01, 1.8385144146980998e+01);
  mol.emplace_back(AtomicNumber(8), 5.8296157034661000e+01, 7.3415854672649999e+01, 1.6451954460233999e+01);
  mol.emplace_back(AtomicNumber(6), 5.3125866728756996e+01, 7.3791910144460999e+01, 1.8925605779834999e+01);
  mol.emplace_back(AtomicNumber(8), 5.1551724979920003e+01, 7.5063695735058005e+01, 2.0779426975044000e+01);
  mol.emplace_back(AtomicNumber(6), 5.2133760584531998e+01, 7.1119837596015003e+01, 1.8358687983134999e+01);
  mol.emplace_back(AtomicNumber(7), 5.8194111831255000e+01, 6.9948207482835002e+01, 1.9076783858955000e+01);
  mol.emplace_back(AtomicNumber(6), 5.9942108371079996e+01, 6.8576266414821006e+01, 1.7340125675064002e+01);
  mol.emplace_back(AtomicNumber(6), 5.8496467989494995e+01, 6.6539141798678997e+01, 1.5985192140951000e+01);
  mol.emplace_back(AtomicNumber(8), 5.6739022819724994e+01, 6.5418534287202007e+01, 1.7083122940559999e+01);
  mol.emplace_back(AtomicNumber(6), 6.2351509007054993e+01, 6.7809037663287000e+01, 1.8772537974725999e+01);
  mol.emplace_back(AtomicNumber(6), 6.2939213789633996e+01, 6.4970669227808997e+01, 1.8594903731759999e+01);
  mol.emplace_back(AtomicNumber(6), 6.2566937769801001e+01, 6.8750121205808995e+01, 2.1609016684215000e+01);
  mol.emplace_back(AtomicNumber(6), 6.5261687030114999e+01, 6.4303595953691996e+01, 2.0257862602079999e+01);
  mol.emplace_back(AtomicNumber(7), 5.9042598800316000e+01, 6.6113953451153989e+01, 1.3600357942833000e+01);
  mol.emplace_back(AtomicNumber(6), 5.7646091294445000e+01, 6.4031475411276006e+01, 1.2305895640368000e+01);
  mol.emplace_back(AtomicNumber(6), 5.9354403588501000e+01, 6.1756245320520001e+01, 1.2181173725094000e+01);
  mol.emplace_back(AtomicNumber(8), 6.1640972035190998e+01, 6.2005689151067997e+01, 1.1574571682625001e+01);
  mol.emplace_back(AtomicNumber(6), 5.6863744734999003e+01, 6.4993345939676999e+01, 9.5960285721420000e+00);
  mol.emplace_back(AtomicNumber(8), 5.9412985094160000e+01, 6.5220113058357001e+01, 8.4792005126430006e+00);
  mol.emplace_back(AtomicNumber(6), 5.5595738596380002e+01, 6.7569042462683996e+01, 9.6735073376909995e+00);
  mol.emplace_back(AtomicNumber(7), 5.8362297444276003e+01, 5.9498022763664999e+01, 1.2596913442674001e+01);
  mol.emplace_back(AtomicNumber(6), 5.9860850153553002e+01, 5.7211454316974994e+01, 1.2545890840971001e+01);
  mol.emplace_back(AtomicNumber(6), 5.8623079630757999e+01, 5.5346294765831999e+01, 1.0705297727685000e+01);
  mol.emplace_back(AtomicNumber(8), 5.6330842006101001e+01, 5.5548495446654996e+01, 1.0478530609005000e+01);
  mol.emplace_back(AtomicNumber(6), 5.9643531664817999e+01, 5.6098405709453999e+01, 1.5202845581505001e+01);
  mol.emplace_back(AtomicNumber(6), 6.1663648747059000e+01, 5.5641092020115998e+01, 1.7120917460339999e+01);
  mol.emplace_back(AtomicNumber(6), 6.3899194592046001e+01, 5.7428772805710004e+01, 1.7064225680669999e+01);
  mol.emplace_back(AtomicNumber(6), 6.0367296718604997e+01, 5.5650540650061004e+01, 1.9721180421204000e+01);
  mol.emplace_back(AtomicNumber(7), 6.0157537133825997e+01, 5.3690894799467998e+01, 9.6848456936250003e+00);
  mol.emplace_back(AtomicNumber(6), 5.8997245376579997e+01, 5.1666998265248999e+01, 8.0785786029750000e+00);
  mol.emplace_back(AtomicNumber(6), 5.9412985094160000e+01, 4.9282164067130999e+01, 9.5998080241199997e+00);
  mol.emplace_back(AtomicNumber(8), 6.1559713817663997e+01, 4.8758709968177996e+01, 1.0319793625929000e+01);
  mol.emplace_back(AtomicNumber(6), 6.0144309051903001e+01, 5.1517709912118001e+01, 5.4688670121660001e+00);
  mol.emplace_back(AtomicNumber(6), 5.9267476193006999e+01, 5.3687115347490000e+01, 3.7076423904180000e+00);
  mol.emplace_back(AtomicNumber(6), 5.9849511797619002e+01, 5.3462237954799001e+01, 9.4108354252199999e-01);
  mol.emplace_back(AtomicNumber(8), 5.8333951554441001e+01, 5.4085847531168994e+01, -6.9163971197399998e-01);
  mol.emplace_back(AtomicNumber(8), 6.2049152848814998e+01, 5.2649655779528999e+01, 5.2534382494200005e-01);
  mol.emplace_back(AtomicNumber(7), 5.7277594726589996e+01, 4.8108644227961996e+01, 1.0174284724776001e+01);
  mol.emplace_back(AtomicNumber(6), 5.7236020754831998e+01, 4.5816406603305005e+01, 1.1703073049876998e+01);
  mol.emplace_back(AtomicNumber(6), 5.5329287231930998e+01, 4.3892665546503004e+01, 1.0659944303949000e+01);
  mol.emplace_back(AtomicNumber(8), 5.3815616714741999e+01, 4.4450134713257995e+01, 8.9289552980249987e+00);
  mol.emplace_back(AtomicNumber(6), 5.6508476249066995e+01, 4.6468362069510000e+01, 1.4484749705684999e+01);
  mol.emplace_back(AtomicNumber(6), 5.8320723472517997e+01, 4.8180453815543999e+01, 1.5852911321720999e+01);
  mol.emplace_back(AtomicNumber(6), 5.3811837262764001e+01, 4.7498262733515006e+01, 1.4560338745245000e+01);
  mol.emplace_back(AtomicNumber(7), 5.5520149556819995e+01, 4.1681686139372999e+01, 1.1776772363448000e+01);
  mol.emplace_back(AtomicNumber(6), 5.3796719454852003e+01, 3.9570862209660000e+01, 1.1300561414220001e+01);
  mol.emplace_back(AtomicNumber(6), 5.2570287287991000e+01, 3.8945362907301003e+01, 1.3825235335523999e+01);
  mol.emplace_back(AtomicNumber(8), 5.3760814661061005e+01, 3.9068195096585995e+01, 1.5798109268039999e+01);
  mol.emplace_back(AtomicNumber(6), 5.5204565316657003e+01, 3.7221932805332997e+01, 1.0404831295434001e+01);
  mol.emplace_back(AtomicNumber(6), 5.6177774200992005e+01, 3.7331536912695000e+01, 7.6722875153399990e+00);
  mol.emplace_back(AtomicNumber(6), 5.4337181087706000e+01, 3.7909793065328998e+01, 5.6276039952420005e+00);
  mol.emplace_back(AtomicNumber(8), 5.2054392092994000e+01, 3.7779401972088003e+01, 5.6408320771649993e+00);
  mol.emplace_back(AtomicNumber(8), 5.5437001613303998e+01, 3.8593873873346993e+01, 3.5980382830559998e+00);
  mol.emplace_back(AtomicNumber(7), 5.0189232541850998e+01, 3.8210259497579997e+01, 1.3772323007832000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8809732569881000e+01, 3.7463817731924998e+01, 1.6051332550565999e+01);
  mol.emplace_back(AtomicNumber(6), 5.0155217474049003e+01, 3.5398347225948001e+01, 1.7481855124238997e+01);
  mol.emplace_back(AtomicNumber(8), 4.9762154468336995e+01, 3.5027960932104001e+01, 1.9760864666973003e+01);
  mol.emplace_back(AtomicNumber(6), 4.6239705224841003e+01, 3.6532182819348002e+01, 1.5027101064528001e+01);
  mol.emplace_back(AtomicNumber(6), 4.5918451806710998e+01, 3.8047743062526003e+01, 1.2668723030256000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8592414081145996e+01, 3.7998610186812002e+01, 1.1476305931197000e+01);
  mol.emplace_back(AtomicNumber(7), 5.1704792785028999e+01, 3.3937589036451001e+01, 1.6174164739850998e+01);
  mol.emplace_back(AtomicNumber(6), 5.3014372895405998e+01, 3.1813537024815002e+01, 1.7404376358690001e+01);
  mol.emplace_back(AtomicNumber(6), 5.5289602986161995e+01, 3.2726274677502005e+01, 1.8867024274176000e+01);
  mol.emplace_back(AtomicNumber(8), 5.6559498850769998e+01, 3.1137015120753002e+01, 2.0042433839333999e+01);
  mol.emplace_back(AtomicNumber(6), 5.3900654384246998e+01, 2.9895465145980001e+01, 1.5461738041998000e+01);
  mol.emplace_back(AtomicNumber(8), 5.4700008477594004e+01, 3.1076543889105000e+01, 1.3165720965362999e+01);
  mol.emplace_back(AtomicNumber(7), 5.5933999548411002e+01, 3.5147013669411002e+01, 1.8572227019891997e+01);
  mol.emplace_back(AtomicNumber(6), 5.8196001557243996e+01, 3.6061641048086997e+01, 1.9966844799774002e+01);
  mol.emplace_back(AtomicNumber(6), 5.7619635130599001e+01, 3.6210929401217996e+01, 2.2752300907559999e+01);
  mol.emplace_back(AtomicNumber(8), 5.5495583118963005e+01, 3.6893120483246996e+01, 2.3510081029149003e+01);
  mol.emplace_back(AtomicNumber(6), 5.8874413187295005e+01, 3.8767728664335003e+01, 1.8987966737472000e+01);
  mol.emplace_back(AtomicNumber(6), 6.0325722746846999e+01, 3.8618440311203997e+01, 1.6544551033695001e+01);
  mol.emplace_back(AtomicNumber(8), 6.1402866560577003e+01, 3.6611551310886000e+01, 1.5979522962983999e+01);
  mol.emplace_back(AtomicNumber(8), 6.0165096037782000e+01, 4.0443915616578003e+01, 1.5057336680352000e+01);
  mol.emplace_back(AtomicNumber(7), 5.9545265913390004e+01, 3.5783851327703999e+01, 2.4286758410628000e+01);
  mol.emplace_back(AtomicNumber(6), 5.9333616602622001e+01, 3.6025736254296000e+01, 2.6996625478854000e+01);
  mol.emplace_back(AtomicNumber(6), 5.9702113170476999e+01, 3.8839538251916998e+01, 2.7693934368794999e+01);
  mol.emplace_back(AtomicNumber(8), 6.0771698080250999e+01, 4.0271950551579003e+01, 2.6193491933529000e+01);
  mol.emplace_back(AtomicNumber(6), 6.1400976834587993e+01, 3.4379784917877004e+01, 2.8336441205054999e+01);
  mol.emplace_back(AtomicNumber(8), 6.3831164456441996e+01, 3.5411575307871004e+01, 2.7431262456323999e+01);
  mol.emplace_back(AtomicNumber(6), 6.1136415196127992e+01, 3.1558424016299998e+01, 2.7646691219070000e+01);
  mol.emplace_back(AtomicNumber(7), 5.8795044695756999e+01, 3.9425353308506999e+01, 2.9971054185539998e+01);
  mol.emplace_back(AtomicNumber(6), 5.9125746743832003e+01, 4.1953806681789004e+01, 3.1023631561413001e+01);
  mol.emplace_back(AtomicNumber(6), 6.1937659015464007e+01, 4.2554739546290996e+01, 3.1325987719653003e+01);
  mol.emplace_back(AtomicNumber(8), 6.2801263792436991e+01, 4.4709027173750997e+01, 3.0961270603776001e+01);
  mol.emplace_back(AtomicNumber(6), 5.7674437184280002e+01, 4.2140889554700003e+01, 3.3569092468595997e+01);
  mol.emplace_back(AtomicNumber(6), 5.4813392036933998e+01, 4.1655229975527000e+01, 3.2960600700138002e+01);
  mol.emplace_back(AtomicNumber(6), 5.8264031692848000e+01, 4.4784616213311004e+01, 3.4691589706062004e+01);
  mol.emplace_back(AtomicNumber(6), 5.3681446169522999e+01, 4.3365431995572003e+01, 3.0927255535973998e+01);
  mol.emplace_back(AtomicNumber(7), 6.3396527478972004e+01, 4.0678241639214001e+01, 3.2030855513550001e+01);
  mol.emplace_back(AtomicNumber(6), 6.6198991120659002e+01, 4.1048627933058000e+01, 3.2255732906241001e+01);
  mol.emplace_back(AtomicNumber(6), 6.7302591098234998e+01, 4.1933019695910005e+01, 2.9780191860651001e+01);
  mol.emplace_back(AtomicNumber(8), 6.9035469830147989e+01, 4.3550625142493999e+01, 2.9714051451035999e+01);
  mol.emplace_back(AtomicNumber(6), 6.7400856849663001e+01, 3.8518284833787000e+01, 3.2970049330083000e+01);
  mol.emplace_back(AtomicNumber(6), 7.0161746519592000e+01, 3.8348209494776995e+01, 3.3773182875407997e+01);
  mol.emplace_back(AtomicNumber(6), 7.0979997872829003e+01, 3.5623224618639000e+01, 3.4170025333098003e+01);
  mol.emplace_back(AtomicNumber(8), 7.1352273892661998e+01, 3.4060421225736000e+01, 3.2493838380855003e+01);
  mol.emplace_back(AtomicNumber(8), 7.1106609514092000e+01, 3.5147013669411002e+01, 3.6496278025556997e+01);
  mol.emplace_back(AtomicNumber(7), 6.6403081527471002e+01, 4.0863434786135997e+01, 2.7707162450718002e+01);
  mol.emplace_back(AtomicNumber(6), 6.7255347948510007e+01, 4.1470036828605004e+01, 2.5137135105677999e+01);
  mol.emplace_back(AtomicNumber(6), 6.6590164400381994e+01, 4.4185573074798000e+01, 2.4415259777879999e+01);
  mol.emplace_back(AtomicNumber(8), 6.8154857519274003e+01, 4.5559403868801006e+01, 2.3305990622336999e+01);
  mol.emplace_back(AtomicNumber(6), 6.6261352078295999e+01, 3.9602987551472999e+01, 2.3158591995195000e+01);
  mol.emplace_back(AtomicNumber(6), 6.7162751375048998e+01, 4.0474151232402001e+01, 2.0543211226419000e+01);
  mol.emplace_back(AtomicNumber(8), 6.9489004067507992e+01, 4.0861545060147002e+01, 2.0174714658564000e+01);
  mol.emplace_back(AtomicNumber(7), 6.5437431547092004e+01, 4.0808632732454996e+01, 1.8746081810880000e+01);
  mol.emplace_back(AtomicNumber(7), 6.4263911707923000e+01, 4.4871543608804998e+01, 2.5038869354250000e+01);
  mol.emplace_back(AtomicNumber(6), 6.3368181589137002e+01, 4.7426453145933003e+01, 2.4524863885241999e+01);
  mol.emplace_back(AtomicNumber(6), 6.5084052787149005e+01, 4.9319958586911000e+01, 2.5859010433475998e+01);
  mol.emplace_back(AtomicNumber(8), 6.5919311674287002e+01, 5.1192677042009997e+01, 2.4742182373976998e+01);
  mol.emplace_back(AtomicNumber(6), 6.0584615207340001e+01, 4.7728809304173005e+01, 2.5254298116996001e+01);
  mol.emplace_back(AtomicNumber(6), 5.9874078235475999e+01, 5.0548280479760997e+01, 2.5212724145238003e+01);
  mol.emplace_back(AtomicNumber(6), 5.8868744009328005e+01, 4.6148998377368997e+01, 2.3578111164753000e+01);
  mol.emplace_back(AtomicNumber(7), 6.5637742501925999e+01, 4.8796504487957996e+01, 2.8249513809560998e+01);
  mol.emplace_back(AtomicNumber(6), 6.7266686304443994e+01, 5.0484029796134998e+01, 2.9736728162904001e+01);
  mol.emplace_back(AtomicNumber(6), 6.9872618443275002e+01, 5.0693789380913998e+01, 2.8548090515822999e+01);
  mol.emplace_back(AtomicNumber(8), 7.1014012940631005e+01, 5.2772487968813998e+01, 2.8646356267251001e+01);
  mol.emplace_back(AtomicNumber(6), 6.7491563697135007e+01, 4.9516490089766997e+01, 3.2450374683108002e+01);
  mol.emplace_back(AtomicNumber(6), 6.4898859640227002e+01, 4.9973803779104998e+01, 3.3822315751121998e+01);
  mol.emplace_back(AtomicNumber(6), 6.5212554154401005e+01, 4.9278384615153001e+01, 3.6585095147040001e+01);
  mol.emplace_back(AtomicNumber(6), 6.2750241190734009e+01, 4.9720580496578997e+01, 3.8025066350658001e+01);
  mol.emplace_back(AtomicNumber(7), 6.3220782961994999e+01, 4.8962800374990003e+01, 4.0716036158994001e+01);
  mol.emplace_back(AtomicNumber(7), 7.0862834861511004e+01, 4.8647216134826998e+01, 2.7535197385718998e+01);
  mol.emplace_back(AtomicNumber(6), 7.3310030017265987e+01, 4.8681231202629000e+01, 2.6229396727320001e+01);
  mol.emplace_back(AtomicNumber(6), 7.3185308101992007e+01, 5.0249703773499000e+01, 2.3831334447279001e+01);
  mol.emplace_back(AtomicNumber(8), 7.5029680667256002e+01, 5.1676446895193997e+01, 2.3200165966952998e+01);
  mol.emplace_back(AtomicNumber(6), 7.4237885477864992e+01, 4.5988371668303998e+01, 2.5636022766774001e+01);
  mol.emplace_back(AtomicNumber(7), 7.1116058144036998e+01, 5.0158996926027001e+01, 2.2425378311463003e+01);
  mol.emplace_back(AtomicNumber(6), 7.0809922533818991e+01, 5.1761484564698996e+01, 2.0159596850651997e+01);
  mol.emplace_back(AtomicNumber(6), 7.0753230754149001e+01, 5.4579066014298000e+01, 2.0885251630427998e+01);
  mol.emplace_back(AtomicNumber(8), 7.1847382101779999e+01, 5.6260922144507994e+01, 1.9619135217798000e+01);
  mol.emplace_back(AtomicNumber(6), 6.8394852719876994e+01, 5.1132205810362002e+01, 1.8729074276978999e+01);
  mol.emplace_back(AtomicNumber(6), 6.8319263680316993e+01, 4.8414779838180003e+01, 1.7780431830501001e+01);
  mol.emplace_back(AtomicNumber(6), 6.5683095925662002e+01, 4.7772273001919999e+01, 1.6818561302100001e+01);
  mol.emplace_back(AtomicNumber(6), 6.5749236335277004e+01, 4.5852311397095995e+01, 1.4677501756563000e+01);
  mol.emplace_back(AtomicNumber(7), 6.5977893179946008e+01, 4.7137325069615997e+01, 1.2171725095149000e+01);
  mol.emplace_back(AtomicNumber(7), 6.9562703381079004e+01, 5.5123307099130002e+01, 2.3039539257887999e+01);
  mol.emplace_back(AtomicNumber(6), 6.9411525301959003e+01, 5.7768923483729999e+01, 2.3895585130904998e+01);
  mol.emplace_back(AtomicNumber(6), 7.2089267028372007e+01, 5.8545600865209003e+01, 2.4696828950241002e+01);
  mol.emplace_back(AtomicNumber(8), 7.2837598520015987e+01, 6.0754690546349998e+01, 2.4294317314583999e+01);
  mol.emplace_back(AtomicNumber(6), 6.7478335615212004e+01, 5.8158207037464003e+01, 2.6089557004133997e+01);
  mol.emplace_back(AtomicNumber(6), 6.4681541151491999e+01, 5.7882307043070000e+01, 2.5169260447491002e+01);
  mol.emplace_back(AtomicNumber(6), 6.7792030129386006e+01, 6.0732013834481997e+01, 2.7423703552368000e+01);
  mol.emplace_back(AtomicNumber(6), 6.2897639817875998e+01, 5.7644201568456005e+01, 2.7499292591928000e+01);
  mol.emplace_back(AtomicNumber(7), 7.3478215630287011e+01, 5.6899649528790000e+01, 2.5913812487156999e+01);
  mol.emplace_back(AtomicNumber(6), 7.6097375851040994e+01, 5.7651760472412001e+01, 2.6673482334734999e+01);
  mol.emplace_back(AtomicNumber(6), 7.7652620339987990e+01, 5.8218678269111997e+01, 2.4284868684639001e+01);
  mol.emplace_back(AtomicNumber(8), 7.9043458667892011e+01, 6.0108404258112003e+01, 2.3963615266508999e+01);
  mol.emplace_back(AtomicNumber(6), 7.7471206645044006e+01, 5.5556054350610999e+01, 2.8090776826485001e+01);
  mol.emplace_back(AtomicNumber(6), 8.0209419603105005e+01, 5.6404541319671999e+01, 2.8689819964998001e+01);
  mol.emplace_back(AtomicNumber(6), 8.1428292866010011e+01, 5.4477020810892000e+01, 3.0415139792954999e+01);
  mol.emplace_back(AtomicNumber(8), 8.0823580549530007e+01, 5.2260372225795003e+01, 3.0057981581033999e+01);
  mol.emplace_back(AtomicNumber(7), 8.2955191465121999e+01, 5.5278264630227994e+01, 3.2219828112450003e+01);
  mol.emplace_back(AtomicNumber(7), 7.7480655274988990e+01, 5.6461233099342003e+01, 2.2546320774759000e+01);
  mol.emplace_back(AtomicNumber(6), 7.8835588809102006e+01, 5.6733353641757994e+01, 2.0112353700927002e+01);
  mol.emplace_back(AtomicNumber(6), 7.8232766218611005e+01, 5.9220233043282001e+01, 1.8834898932363000e+01);
  mol.emplace_back(AtomicNumber(8), 7.9859820295139997e+01, 6.0539261783604005e+01, 1.7727519502808999e+01);
  mol.emplace_back(AtomicNumber(6), 7.8230876492622002e+01, 5.4386313963420001e+01, 1.8538211952090002e+01);
  mol.emplace_back(AtomicNumber(6), 8.0551460007114002e+01, 5.3964905067872998e+01, 1.6871473629792000e+01);
  mol.emplace_back(AtomicNumber(8), 8.2516775035673987e+01, 5.3407435901117999e+01, 1.8026096209071000e+01);
  mol.emplace_back(AtomicNumber(8), 8.0181073713269996e+01, 5.4446785195068003e+01, 1.4603802442992000e+01);
  mol.emplace_back(AtomicNumber(7), 7.5810137500712997e+01, 5.9998800150750000e+01, 1.8874583178131999e+01);
  mol.emplace_back(AtomicNumber(6), 7.5226212170112007e+01, 6.2349619281065998e+01, 1.7447840056437002e+01);
  mol.emplace_back(AtomicNumber(6), 7.5281014223793008e+01, 6.4762799369019007e+01, 1.8887811260054999e+01);
  mol.emplace_back(AtomicNumber(8), 7.5898954622196001e+01, 6.6750791109447007e+01, 1.7659489367205001e+01);
  mol.emplace_back(AtomicNumber(6), 7.2971769065235009e+01, 6.1984902165189006e+01, 1.5722520228480001e+01);
  mol.emplace_back(AtomicNumber(6), 7.0335601310580003e+01, 6.2024586410958008e+01, 1.6680611304903000e+01);
  mol.emplace_back(AtomicNumber(6), 6.8693429426139005e+01, 6.3519359668256996e+01, 1.4811672301782000e+01);
  mol.emplace_back(AtomicNumber(6), 6.8638627372458004e+01, 6.2255132981616008e+01, 1.2239755230753001e+01);
  mol.emplace_back(AtomicNumber(7), 6.7591719174551997e+01, 6.4146748696605002e+01, 1.0372705953620999e+01);
  mol.emplace_back(AtomicNumber(7), 7.4937084093795008e+01, 6.4883741832314996e+01, 2.1325557785865001e+01);
  mol.emplace_back(AtomicNumber(6), 7.4976768339564003e+01, 6.7174089730982999e+01, 2.2812772139207997e+01);
  mol.emplace_back(AtomicNumber(6), 7.6864604602575000e+01, 6.7136295211203006e+01, 2.4944383054799999e+01);
  mol.emplace_back(AtomicNumber(8), 7.7127276515046006e+01, 6.9027910926191993e+01, 2.6287978232979000e+01);
  mol.emplace_back(AtomicNumber(6), 7.2357608118809992e+01, 6.7678646570045998e+01, 2.3995740608322002e+01);
  mol.emplace_back(AtomicNumber(6), 7.0214658847283999e+01, 6.8001789714165000e+01, 2.2087117359432000e+01);
  mol.emplace_back(AtomicNumber(6), 7.0282688982888004e+01, 7.0602052675029000e+01, 2.0849346836637000e+01);
  mol.emplace_back(AtomicNumber(8), 7.0900629381290997e+01, 7.2489888938039996e+01, 2.2005859141904999e+01);
  mol.emplace_back(AtomicNumber(8), 6.9657189680529001e+01, 7.0524573909479997e+01, 1.8560888663958000e+01);
  mol.emplace_back(AtomicNumber(7), 7.8077808687512999e+01, 6.4993345939676999e+01, 2.5382799484248000e+01);
  mol.emplace_back(AtomicNumber(6), 8.0020447004204996e+01, 6.4759019917041002e+01, 2.7270635747259000e+01);
  mol.emplace_back(AtomicNumber(6), 7.9272115512561001e+01, 6.4394302801164002e+01, 2.9937039117737999e+01);
  mol.emplace_back(AtomicNumber(8), 8.0935074382880998e+01, 6.4250683625999997e+01, 3.1632123329871000e+01);
  mol.emplace_back(AtomicNumber(7), 7.6802243644938002e+01, 6.4091946642924000e+01, 3.0447265134767996e+01);
  mol.emplace_back(AtomicNumber(6), 7.6016117633514000e+01, 6.3714001445124005e+01, 3.3087212341400999e+01);
  mol.emplace_back(AtomicNumber(6), 7.6437526529061003e+01, 6.0996575472941998e+01, 3.3911132872605002e+01);
  mol.emplace_back(AtomicNumber(8), 7.5468097096704000e+01, 5.9216453591303996e+01, 3.2720605499535004e+01);
  mol.emplace_back(AtomicNumber(6), 7.3119167692376990e+01, 6.4450994580834006e+01, 3.3249728776455001e+01);
  mol.emplace_back(AtomicNumber(6), 7.2699648522818990e+01, 6.7172200004993996e+01, 3.2210379482505004e+01);
  mol.emplace_back(AtomicNumber(6), 7.2085487576394002e+01, 6.4122182258747998e+01, 3.5955816392702999e+01);
  mol.emplace_back(AtomicNumber(6), 6.9840493101462002e+01, 6.7550145202793999e+01, 3.1520629496519998e+01);
  mol.emplace_back(AtomicNumber(7), 7.7835923760921006e+01, 6.0631858357064999e+01, 3.5963375296658995e+01);
  mol.emplace_back(AtomicNumber(6), 7.8349929229929003e+01, 5.8110963887739004e+01, 3.7027291028466003e+01);
  mol.emplace_back(AtomicNumber(6), 7.5906513526151997e+01, 5.6740912545713996e+01, 3.7639562248901996e+01);
  mol.emplace_back(AtomicNumber(8), 7.4198201232096011e+01, 5.7942778274717995e+01, 3.8779067020269004e+01);
  mol.emplace_back(AtomicNumber(6), 7.9736988105855005e+01, 5.8849846749438001e+01, 3.9519839607957003e+01);
  mol.emplace_back(AtomicNumber(6), 8.1076803832056001e+01, 6.1253578207446004e+01, 3.8839538251916998e+01);
  mol.emplace_back(AtomicNumber(6), 7.9032120311958010e+01, 6.2716226122932007e+01, 3.7441141020056996e+01);
  mol.emplace_back(AtomicNumber(7), 7.5700533393350995e+01, 5.4344739991661996e+01, 3.7051857466323000e+01);
  mol.emplace_back(AtomicNumber(6), 7.3353493715013002e+01, 5.2950122211779998e+01, 3.7584760195221001e+01);
  mol.emplace_back(AtomicNumber(6), 7.2605162223368993e+01, 5.3003034539471997e+01, 4.0328642331249000e+01);
  mol.emplace_back(AtomicNumber(8), 7.0322373228657000e+01, 5.2980357827604003e+01, 4.1014612865255998e+01);
  mol.emplace_back(AtomicNumber(6), 7.3869388910010002e+01, 5.0321513361081003e+01, 3.6518954737424998e+01);
  mol.emplace_back(AtomicNumber(6), 7.5743997091097995e+01, 5.0841188008056001e+01, 3.4389233547822002e+01);
  mol.emplace_back(AtomicNumber(6), 7.7544905958614990e+01, 5.2740362627000998e+01, 3.5676136946330999e+01);
  mol.emplace_back(AtomicNumber(7), 7.4406071090886002e+01, 5.3082403031010003e+01, 4.2027505995359995e+01);
  mol.emplace_back(AtomicNumber(6), 7.3818366308307006e+01, 5.3031380429306999e+01, 4.4777057309355001e+01);
  mol.emplace_back(AtomicNumber(6), 7.2499337567985009e+01, 5.5435111887315003e+01, 4.5653890168250996e+01);
  mol.emplace_back(AtomicNumber(8), 7.1212434169475998e+01, 5.5539046816709998e+01, 4.7660779168569000e+01);
  mol.emplace_back(AtomicNumber(6), 7.6231546396260001e+01, 5.2330292087388003e+01, 4.6237815498852001e+01);
  mol.emplace_back(AtomicNumber(6), 7.6645396387850994e+01, 5.4017817395565004e+01, 4.8518714767574998e+01);
  mol.emplace_back(AtomicNumber(8), 7.6942083368124003e+01, 5.6330842006101001e+01, 4.8104864775983998e+01);
  mol.emplace_back(AtomicNumber(8), 7.6626499127960997e+01, 5.3082403031010003e+01, 5.0720245544759997e+01);
  mol.emplace_back(AtomicNumber(7), 7.2601382771390988e+01, 5.7396647463897004e+01, 4.4108094309249005e+01);
  mol.emplace_back(AtomicNumber(6), 7.1314479372882005e+01, 5.9785261113993002e+01, 4.4809182651168001e+01);
  mol.emplace_back(AtomicNumber(6), 6.8661304084326005e+01, 5.9983682342838001e+01, 4.3628103908043002e+01);
  mol.emplace_back(AtomicNumber(8), 6.7225112332685995e+01, 6.1639082309202003e+01, 4.4376435399687004e+01);
  mol.emplace_back(AtomicNumber(6), 7.2807362904192004e+01, 6.2085057642605996e+01, 4.3807627876997998e+01);
  mol.emplace_back(AtomicNumber(6), 7.5435971754891000e+01, 6.2085057642605996e+01, 4.5051067577760001e+01);
  mol.emplace_back(AtomicNumber(6), 7.7025231311639999e+01, 6.4318713761604002e+01, 4.4208249786665995e+01);
  mol.emplace_back(AtomicNumber(8), 7.9321248388275009e+01, 6.4265801433912003e+01, 4.4642886764135994e+01);
  mol.emplace_back(AtomicNumber(7), 7.5853601198459998e+01, 6.6153637696922999e+01, 4.3038509399474997e+01);
  mol.emplace_back(AtomicNumber(7), 6.8030135603999994e+01, 5.8316944020539999e+01, 4.1899004628108003e+01);
  mol.emplace_back(AtomicNumber(6), 6.5645301405881995e+01, 5.8345289910375001e+01, 4.0578086161796996e+01);
  mol.emplace_back(AtomicNumber(6), 6.3474006244521000e+01, 5.7048937881920999e+01, 4.1916012162009004e+01);
  mol.emplace_back(AtomicNumber(8), 6.3456998710619999e+01, 5.4819061214900998e+01, 4.2516945026510996e+01);
  mol.emplace_back(AtomicNumber(6), 6.5906083592363998e+01, 5.7139644729392998e+01, 3.7919241695273996e+01);
  mol.emplace_back(AtomicNumber(6), 6.8052812315867996e+01, 5.8316944020539999e+01, 3.6322423234569001e+01);
  mol.emplace_back(AtomicNumber(6), 6.8186982861087003e+01, 5.7058386511865997e+01, 3.3778852053374997e+01);
  mol.emplace_back(AtomicNumber(8), 6.6231116462472002e+01, 5.6128641325278004e+01, 3.2868004126677000e+01);
  mol.emplace_back(AtomicNumber(7), 7.0350719118492009e+01, 5.6929885144614005e+01, 3.2565647968436998e+01);
  mol.emplace_back(AtomicNumber(7), 6.1374520670742001e+01, 5.8424658401913000e+01, 4.2082308049040996e+01);
  mol.emplace_back(AtomicNumber(6), 5.8959450856799997e+01, 5.7313499520381001e+01, 4.3047958029420002e+01);
  mol.emplace_back(AtomicNumber(6), 5.7088622127690002e+01, 5.7653650198401003e+01, 4.0912567661849998e+01);
  mol.emplace_back(AtomicNumber(8), 5.6650205698242004e+01, 5.9953446727013997e+01, 4.0192582060040998e+01);
  mol.emplace_back(AtomicNumber(6), 5.8292377582683002e+01, 5.8451114565758999e+01, 4.5576411402702000e+01);
  mol.emplace_back(AtomicNumber(6), 5.5580620788467996e+01, 5.8196001557243996e+01, 4.6483479877421999e+01);
  mol.emplace_back(AtomicNumber(6), 5.5314169424018999e+01, 5.9174879619545997e+01, 4.9163111329823998e+01);
  mol.emplace_back(AtomicNumber(7), 5.2676111943374998e+01, 5.9180548797512998e+01, 4.9970024327127000e+01);
  mol.emplace_back(AtomicNumber(6), 5.1272045533548003e+01, 6.1270585741347006e+01, 5.0217578431686000e+01);
  mol.emplace_back(AtomicNumber(7), 5.2213129076069997e+01, 6.3600617885783997e+01, 5.0004039394928995e+01);
  mol.emplace_back(AtomicNumber(7), 4.8773827776089995e+01, 6.1036259718711001e+01, 5.0516155137947997e+01);
  mol.emplace_back(AtomicNumber(7), 5.6113523517365998e+01, 5.5625974212204000e+01, 3.9786290972406000e+01);
  mol.emplace_back(AtomicNumber(6), 5.4352298895617999e+01, 5.5884866672697001e+01, 3.7616885537034001e+01);
  mol.emplace_back(AtomicNumber(6), 5.1648101005358996e+01, 5.5401096819513000e+01, 3.8482380039996002e+01);
  mol.emplace_back(AtomicNumber(8), 5.1213464027888996e+01, 5.3566172884194003e+01, 3.9867549189933001e+01);
  mol.emplace_back(AtomicNumber(6), 5.5087402305338998e+01, 5.4150098214795001e+01, 3.5441810923694995e+01);
  mol.emplace_back(AtomicNumber(6), 5.7477905681423998e+01, 5.4635757793967997e+01, 3.3977273282220004e+01);
  mol.emplace_back(AtomicNumber(6), 5.8086397449882000e+01, 5.2332181813377005e+01, 3.2355888383657998e+01);
  mol.emplace_back(AtomicNumber(6), 5.7079173497744996e+01, 5.7009253636151996e+01, 3.2369116465581001e+01);
  mol.emplace_back(AtomicNumber(7), 4.9956796245203996e+01, 5.7130196099448000e+01, 3.7802078683956005e+01);
  mol.emplace_back(AtomicNumber(6), 4.7307400408625995e+01, 5.7013033088130001e+01, 3.8552299901589002e+01);
  mol.emplace_back(AtomicNumber(6), 4.5544286060889000e+01, 5.6973348842360998e+01, 3.6275180084844003e+01);
  mol.emplace_back(AtomicNumber(8), 4.5723810029844003e+01, 5.8483239907571999e+01, 3.4557419160842997e+01);
  mol.emplace_back(AtomicNumber(6), 4.6560958642971002e+01, 5.9386528930314000e+01, 4.0224707401854005e+01);
  mol.emplace_back(AtomicNumber(6), 4.8463912713894004e+01, 5.9847622071630006e+01, 4.2369546399369000e+01);
  mol.emplace_back(AtomicNumber(6), 4.3805738151009002e+01, 5.9165430989600999e+01, 4.1241379983936000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8713356544442000e+01, 5.7515700201203998e+01, 4.4136440199084007e+01);
  mol.emplace_back(AtomicNumber(7), 4.3730149111448995e+01, 5.5155432440943002e+01, 3.6360217754349001e+01);
  mol.emplace_back(AtomicNumber(6), 4.1812077232614001e+01, 5.4919216692318003e+01, 3.4360887657987000e+01);
  mol.emplace_back(AtomicNumber(6), 3.9372440980815000e+01, 5.4100965339081000e+01, 3.5723380096055998e+01);
  mol.emplace_back(AtomicNumber(8), 3.9345984816969001e+01, 5.2409660578926001e+01, 3.7320198556761000e+01);
  mol.emplace_back(AtomicNumber(6), 4.2507496396565998e+01, 5.3020042073372998e+01, 3.2331321945801001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0528953286082995e+01, 5.2664773587440997e+01, 3.0284748699714001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0298406715424996e+01, 5.4448674921056998e+01, 2.8355338464945003e+01);
  mol.emplace_back(AtomicNumber(6), 3.9000164960982005e+01, 5.0521824315914998e+01, 3.0335771301417001e+01);
  mol.emplace_back(AtomicNumber(6), 3.8491828669941000e+01, 5.4136870132871998e+01, 2.6458053571988998e+01);
  mol.emplace_back(AtomicNumber(6), 3.7184138285552997e+01, 5.0151438022071005e+01, 2.8442265860439001e+01);
  mol.emplace_back(AtomicNumber(6), 3.7025401302477000e+01, 5.1901324287884997e+01, 2.6495848091769002e+01);
  mol.emplace_back(AtomicNumber(7), 3.7435471842089996e+01, 5.5516370104841997e+01, 3.5107329423642000e+01);
  mol.emplace_back(AtomicNumber(6), 3.4852216415127003e+01, 5.5072284497426999e+01, 3.6061641048086997e+01);
  mol.emplace_back(AtomicNumber(6), 3.4871113675017000e+01, 5.4690559847648998e+01, 3.8911347839499001e+01);
  mol.emplace_back(AtomicNumber(8), 3.3750506163539995e+01, 5.2900989336065997e+01, 3.9926130695592001e+01);
  mol.emplace_back(AtomicNumber(6), 3.3758065067495998e+01, 5.2868863994252997e+01, 3.4668912994194002e+01);
  mol.emplace_back(AtomicNumber(7), 3.6229826661108000e+01, 5.6328952280111999e+01, 4.0143449184326997e+01);
  mol.emplace_back(AtomicNumber(6), 3.6658794460610999e+01, 5.6491468715165993e+01, 4.2811742280795002e+01);
  mol.emplace_back(AtomicNumber(6), 3.7951367037086996e+01, 5.4289937937981001e+01, 4.4070299789469004e+01);
  mol.emplace_back(AtomicNumber(8), 3.7777512246099000e+01, 5.4015927669576001e+01, 4.6413560015828999e+01);
  mol.emplace_back(AtomicNumber(7), 3.9308190297189000e+01, 5.2781936598759003e+01, 4.2666233379642001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0723595062950004e+01, 5.0637097601244001e+01, 4.3715031303536996e+01);
  mol.emplace_back(AtomicNumber(6), 4.3550625142493999e+01, 5.1187007864042997e+01, 4.3299291585957000e+01);
  mol.emplace_back(AtomicNumber(8), 4.4187462800786996e+01, 5.2207459898102996e+01, 4.1328307379430001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0141559458338001e+01, 4.8223917513290999e+01, 4.2312854619698996e+01);
  mol.emplace_back(AtomicNumber(6), 3.7344764994618004e+01, 4.7388658626153003e+01, 4.2433797082995000e+01);
  mol.emplace_back(AtomicNumber(6), 3.7102880068026003e+01, 4.5136105247265000e+01, 4.0687690269158999e+01);
  mol.emplace_back(AtomicNumber(6), 3.5509841059298999e+01, 4.5771053179569002e+01, 3.8386004014556995e+01);
  mol.emplace_back(AtomicNumber(7), 3.2956821248160004e+01, 4.6591194258794999e+01, 3.9357323172903001e+01);
  mol.emplace_back(AtomicNumber(7), 4.5126656617319995e+01, 5.0506706508002999e+01, 4.5071854563639000e+01);
  mol.emplace_back(AtomicNumber(6), 4.7902664095161001e+01, 5.0780716776407999e+01, 4.4678791557926999e+01);
  mol.emplace_back(AtomicNumber(6), 4.8647216134826998e+01, 4.8350529154553996e+01, 4.3316299119858002e+01);
  mol.emplace_back(AtomicNumber(8), 4.7857310671424997e+01, 4.6277499744621004e+01, 4.4178014170841998e+01);
  mol.emplace_back(AtomicNumber(6), 4.9265156533229998e+01, 5.1069844852724998e+01, 4.7167560685440002e+01);
  mol.emplace_back(AtomicNumber(6), 5.2067620174917003e+01, 5.1695344155084001e+01, 4.6666783298355000e+01);
  mol.emplace_back(AtomicNumber(6), 5.3407435901117999e+01, 5.2111083872664004e+01, 4.9170670233780001e+01);
  mol.emplace_back(AtomicNumber(8), 5.5159211892921000e+01, 5.0720245544759997e+01, 4.9765933920315000e+01);
  mol.emplace_back(AtomicNumber(7), 5.2490918796453002e+01, 5.4017817395565004e+01, 5.0529383219871001e+01);
  mol.emplace_back(AtomicNumber(7), 5.0011598298884998e+01, 4.8545170931420998e+01, 4.1258387517836994e+01);
  mol.emplace_back(AtomicNumber(6), 5.0693789380913998e+01, 4.6337970976268998e+01, 3.9706922480868002e+01);
  mol.emplace_back(AtomicNumber(6), 5.2900989336065997e+01, 4.4939573744408996e+01, 4.0899339579927002e+01);
  mol.emplace_back(AtomicNumber(8), 5.4620639986055998e+01, 4.6192462075115998e+01, 4.1759164904921995e+01);
  mol.emplace_back(AtomicNumber(6), 5.1103859920527000e+01, 4.7228031917088003e+01, 3.6983827330719002e+01);
  mol.emplace_back(AtomicNumber(6), 4.9002484620758999e+01, 4.8838078459716002e+01, 3.5827315025451000e+01);
  mol.emplace_back(AtomicNumber(6), 4.9516490089766997e+01, 4.9289722971086995e+01, 3.3015402753819004e+01);
  mol.emplace_back(AtomicNumber(6), 4.6443795631653003e+01, 4.7602197662910001e+01, 3.6054082144131002e+01);
  mol.emplace_back(AtomicNumber(7), 5.2802723584638002e+01, 4.2420569001072003e+01, 4.0908788209872000e+01);
  mol.emplace_back(AtomicNumber(6), 5.4830399570834999e+01, 4.0925795743773001e+01, 4.2118212842832001e+01);
  mol.emplace_back(AtomicNumber(6), 5.6582175562638000e+01, 3.9884556723834002e+01, 4.0137780006359996e+01);
  mol.emplace_back(AtomicNumber(8), 5.5690224895829999e+01, 3.9073864274552996e+01, 3.8153567717910001e+01);
  mol.emplace_back(AtomicNumber(6), 5.3569952336172001e+01, 3.8814971814060002e+01, 4.3588419662273999e+01);
  mol.emplace_back(AtomicNumber(6), 5.5268816000283003e+01, 3.6766508841983999e+01, 4.4795954569244998e+01);
  mol.emplace_back(AtomicNumber(6), 5.4276709856058005e+01, 3.5993610912483000e+01, 4.7367871640273997e+01);
  mol.emplace_back(AtomicNumber(8), 5.5064725593470996e+01, 3.4264511632548000e+01, 4.8652885312793998e+01);
  mol.emplace_back(AtomicNumber(8), 5.2490918796453002e+01, 3.7495943073737998e+01, 4.7936679162963003e+01);
  mol.emplace_back(AtomicNumber(7), 5.9021811814437001e+01, 3.9854321108009998e+01, 4.0551629997950997e+01);
  mol.emplace_back(AtomicNumber(6), 6.0966339857118001e+01, 3.9060636192630000e+01, 3.8765838938346000e+01);
  mol.emplace_back(AtomicNumber(6), 6.0713116574592000e+01, 3.6592654050996003e+01, 3.7322088282750002e+01);
  mol.emplace_back(AtomicNumber(8), 6.1503022037994000e+01, 3.6503836929513000e+01, 3.5069534903861999e+01);
  mol.emplace_back(AtomicNumber(6), 6.3566602817981995e+01, 3.9147563588124001e+01, 4.0141559458338001e+01);
  mol.emplace_back(AtomicNumber(6), 6.4579495948086006e+01, 4.1817746410581002e+01, 4.0353208769105997e+01);
  mol.emplace_back(AtomicNumber(8), 6.6616620564228000e+01, 4.2182463526458001e+01, 4.1494603266461993e+01);
  mol.emplace_back(AtomicNumber(8), 6.3388968575015994e+01, 4.3626214182053999e+01, 3.9463147828286999e+01);
  mol.emplace_back(AtomicNumber(7), 5.9898644673332996e+01, 3.4602772584579000e+01, 3.8561748531534001e+01);
  mol.emplace_back(AtomicNumber(6), 5.9654870020752000e+01, 3.2053532225417996e+01, 3.7463817731924998e+01);
  mol.emplace_back(AtomicNumber(6), 5.7296491986480000e+01, 3.1554644564322000e+01, 3.6001169816438995e+01);
  mol.emplace_back(AtomicNumber(8), 5.7065945415822000e+01, 2.9587439809772999e+01, 3.4706707513974003e+01);
  mol.emplace_back(AtomicNumber(7), 5.5444560517260001e+01, 3.3247839050466006e+01, 3.6048412966164001e+01);
  mol.emplace_back(AtomicNumber(6), 5.3116418098811998e+01, 3.2954931522171002e+01, 3.4536632174963998e+01);
  mol.emplace_back(AtomicNumber(6), 5.3620974937874998e+01, 3.4013178076010995e+01, 3.1911802776243000e+01);
  mol.emplace_back(AtomicNumber(8), 5.5418104353414002e+01, 3.5500392429354001e+01, 3.1539526756410002e+01);
  mol.emplace_back(AtomicNumber(6), 5.0882761979813999e+01, 3.4376005465898999e+01, 3.5700703384187996e+01);
  mol.emplace_back(AtomicNumber(6), 5.0306395553168997e+01, 3.3635232878210999e+01, 3.8459703328128001e+01);
  mol.emplace_back(AtomicNumber(6), 4.9151772973890004e+01, 3.0934814439930001e+01, 3.8323643056920005e+01);
  mol.emplace_back(AtomicNumber(7), 5.0975358553275001e+01, 2.9330437075269000e+01, 3.9574641661637997e+01);
  mol.emplace_back(AtomicNumber(6), 5.2162106474367000e+01, 2.7255517939347001e+01, 3.9032290302795005e+01);
  mol.emplace_back(AtomicNumber(7), 5.1927780451730996e+01, 2.5951607006937000e+01, 3.6919576647092995e+01);
  mol.emplace_back(AtomicNumber(7), 5.3893095480290995e+01, 2.6393802888363002e+01, 4.0723595062950004e+01);
  mol.emplace_back(AtomicNumber(7), 5.1986361957390002e+01, 3.3427363019421001e+01, 3.0148688428506002e+01);
  mol.emplace_back(AtomicNumber(6), 5.2107304420686006e+01, 3.4377895191888001e+01, 2.7520079577807000e+01);
  mol.emplace_back(AtomicNumber(6), 5.0043723640697998e+01, 3.6433917067920000e+01, 2.7272525473247999e+01);
  mol.emplace_back(AtomicNumber(8), 4.8393992852301004e+01, 3.6638007474732007e+01, 2.8888241193843001e+01);
  mol.emplace_back(AtomicNumber(6), 5.1587629773711001e+01, 3.2229276742395001e+01, 2.5573661809137000e+01);
  mol.emplace_back(AtomicNumber(8), 4.8991146264824998e+01, 3.1390238403279000e+01, 2.6291757684957002e+01);
  mol.emplace_back(AtomicNumber(6), 5.3358303025403998e+01, 2.9978613089496001e+01, 2.5620904958861999e+01);
  mol.emplace_back(AtomicNumber(7), 5.0238365417564999e+01, 3.7913572517306996e+01, 2.5280754280842000e+01);
  mol.emplace_back(AtomicNumber(6), 4.8365646962466002e+01, 3.9890225901801003e+01, 2.4702498128207999e+01);
  mol.emplace_back(AtomicNumber(6), 4.5808847699349002e+01, 3.8618440311203997e+01, 2.4296207040572998e+01);
  mol.emplace_back(AtomicNumber(8), 4.3962585408095997e+01, 3.9591649195538999e+01, 2.5188157707381002e+01);
  mol.emplace_back(AtomicNumber(6), 4.9291612697075998e+01, 4.1362322447232003e+01, 2.2361127627837000e+01);
  mol.emplace_back(AtomicNumber(6), 5.1827624974313999e+01, 4.2738042967223997e+01, 2.2491518721077998e+01);
  mol.emplace_back(AtomicNumber(6), 5.2379424963101997e+01, 4.4108094309249005e+01, 1.9989521511642000e+01);
  mol.emplace_back(AtomicNumber(6), 5.1740697578819997e+01, 4.4826190185069002e+01, 2.4481400187495002e+01);
  mol.emplace_back(AtomicNumber(7), 4.5806957973359999e+01, 3.6345099946437003e+01, 2.3141584461294002e+01);
  mol.emplace_back(AtomicNumber(6), 4.3320078571836000e+01, 3.5116778053586998e+01, 2.2723955017725000e+01);
  mol.emplace_back(AtomicNumber(6), 4.2006719009480996e+01, 3.4476160943315996e+01, 2.5180598803424999e+01);
  mol.emplace_back(AtomicNumber(8), 3.9614325907407000e+01, 3.4493168477216997e+01, 2.5312879622655000e+01);
  mol.emplace_back(AtomicNumber(6), 4.3575191580351003e+01, 3.2741392485414003e+01, 2.1078003681306001e+01);
  mol.emplace_back(AtomicNumber(8), 4.5190907300946002e+01, 3.0982057589655000e+01, 2.2213729000695000e+01);
  mol.emplace_back(AtomicNumber(7), 4.3458028569032997e+01, 3.3973493830242006e+01, 2.7147803557974001e+01);
  mol.emplace_back(AtomicNumber(6), 4.2363877221401999e+01, 3.3330986993982002e+01, 2.9655469945377000e+01);
  mol.emplace_back(AtomicNumber(6), 4.0553519723939999e+01, 3.5407795855892999e+01, 3.0543641160206999e+01);
  mol.emplace_back(AtomicNumber(8), 3.8733713596533001e+01, 3.4971269152433997e+01, 3.1936369214099997e+01);
  mol.emplace_back(AtomicNumber(6), 4.4334861427928999e+01, 3.2750841115359002e+01, 3.1635902781848998e+01);
  mol.emplace_back(AtomicNumber(6), 4.5701133317976002e+01, 3.0265851439823997e+01, 3.1405356211190998e+01);
  mol.emplace_back(AtomicNumber(8), 4.7815736699666999e+01, 3.0035304869166001e+01, 3.2412580163328002e+01);
  mol.emplace_back(AtomicNumber(8), 4.4544621012707999e+01, 2.8548090515822999e+01, 3.0188372674274998e+01);
  mol.emplace_back(AtomicNumber(7), 4.1282953955693998e+01, 3.7707592384506000e+01, 3.0056091855045000e+01);
  mol.emplace_back(AtomicNumber(6), 3.9833534122130999e+01, 3.9965814941361003e+01, 3.0709937047239002e+01);
  mol.emplace_back(AtomicNumber(6), 3.8062860870438001e+01, 4.0799184102509997e+01, 2.8627459007360997e+01);
  mol.emplace_back(AtomicNumber(8), 3.6847767059511000e+01, 4.2792845020904998e+01, 2.8952491877469001e+01);
  mol.emplace_back(AtomicNumber(6), 4.1734598467064998e+01, 4.2053962159206002e+01, 3.1333546623608999e+01);
  mol.emplace_back(AtomicNumber(6), 4.3359762817605002e+01, 4.1481375184538997e+01, 3.3608776714365000e+01);
  mol.emplace_back(AtomicNumber(6), 4.5867429205007994e+01, 4.0712256707016003e+01, 3.3342325349915996e+01);
  mol.emplace_back(AtomicNumber(6), 4.2399782015193004e+01, 4.1870658738273001e+01, 3.6027625980285002e+01);
  mol.emplace_back(AtomicNumber(6), 4.7341415476427997e+01, 4.0222817675865002e+01, 3.5481495169463997e+01);
  mol.emplace_back(AtomicNumber(6), 4.3849201848756003e+01, 4.1398227241023001e+01, 3.8157347169887998e+01);
  mol.emplace_back(AtomicNumber(6), 4.6330412072312996e+01, 4.0572416983829996e+01, 3.7851211559669999e+01);
  mol.emplace_back(AtomicNumber(8), 4.7711801770272004e+01, 4.0254943017678002e+01, 4.0045183432898995e+01);
  mol.emplace_back(AtomicNumber(7), 3.7781291698076998e+01, 3.9465037554276002e+01, 2.6548760419461001e+01);
  mol.emplace_back(AtomicNumber(6), 3.6027625980285002e+01, 4.0349429317127999e+01, 2.4564548131011001e+01);
  mol.emplace_back(AtomicNumber(6), 3.6740052678137999e+01, 4.2981817619805000e+01, 2.3640472122390001e+01);
  mol.emplace_back(AtomicNumber(8), 3.5094101341719004e+01, 4.4616430600290002e+01, 2.3222842678820999e+01);
  mol.emplace_back(AtomicNumber(6), 3.3232721242554000e+01, 4.0217148497898002e+01, 2.5437601537929002e+01);
  mol.emplace_back(AtomicNumber(6), 3.1324097993664001e+01, 4.0171795074161999e+01, 2.3271975554535000e+01);
  mol.emplace_back(AtomicNumber(8), 2.9177369270160000e+01, 4.1231931353991001e+01, 2.3391028291842002e+01);
  mol.emplace_back(AtomicNumber(7), 3.1981722637836000e+01, 3.8901899209553996e+01, 2.1195166692623999e+01);
  mol.emplace_back(AtomicNumber(7), 3.9149453314112996e+01, 4.3395667611396000e+01, 2.3168040625139998e+01);
  mol.emplace_back(AtomicNumber(6), 4.0031955350975998e+01, 4.5850421671107000e+01, 2.2090896811409998e+01);
  mol.emplace_back(AtomicNumber(6), 3.9892115627789998e+01, 4.5563183320779004e+01, 1.9224182486097000e+01);
  mol.emplace_back(AtomicNumber(8), 4.1273505325748999e+01, 4.3837863492822002e+01, 1.8303885929454001e+01);
  mol.emplace_back(AtomicNumber(6), 4.2802293650849997e+01, 4.6328522346323993e+01, 2.3001744738108002e+01);
  mol.emplace_back(AtomicNumber(6), 4.2824970362717998e+01, 4.6901109320990997e+01, 2.5887356323311000e+01);
  mol.emplace_back(AtomicNumber(6), 4.4174234718864000e+01, 4.8462022987905002e+01, 2.1559883808501002e+01);
  mol.emplace_back(AtomicNumber(6), 4.5585860032647005e+01, 4.7207244931209004e+01, 2.6824660413855000e+01);
  mol.emplace_back(AtomicNumber(7), 3.8344430042798997e+01, 4.7006933976375002e+01, 1.7965624977422998e+01);
  mol.emplace_back(AtomicNumber(6), 3.7947587585108998e+01, 4.6814181925496996e+01, 1.5180168869636999e+01);
  mol.emplace_back(AtomicNumber(6), 3.9347874542957996e+01, 4.8970359278946006e+01, 1.3855470951348000e+01);
  mol.emplace_back(AtomicNumber(8), 4.0294627263446998e+01, 5.0701348284869994e+01, 1.5132925719911999e+01);
  mol.emplace_back(AtomicNumber(6), 3.5147013669411002e+01, 4.6744262063904003e+01, 1.4601912717003001e+01);
  mol.emplace_back(AtomicNumber(6), 3.3673027397991000e+01, 4.4283838826226003e+01, 1.4928835313100000e+01);
  mol.emplace_back(AtomicNumber(6), 3.1197486352401000e+01, 4.4463362795180998e+01, 1.3447290137724000e+01);
  mol.emplace_back(AtomicNumber(8), 2.9188707626094001e+01, 4.3425903227219997e+01, 1.4046333276237000e+01);
  mol.emplace_back(AtomicNumber(7), 3.1254178132071004e+01, 4.5907113450776997e+01, 1.1355363467901000e+01);
  mol.emplace_back(AtomicNumber(7), 3.9540626593836002e+01, 4.8872093527517997e+01, 1.1349694289934000e+01);
  mol.emplace_back(AtomicNumber(6), 4.0923906017783999e+01, 5.0733473626683001e+01, 9.9021641823600000e+00);
  mol.emplace_back(AtomicNumber(6), 3.9924240969602998e+01, 5.3365861929359994e+01, 1.0533332662686000e+01);
  mol.emplace_back(AtomicNumber(8), 3.7715151288461996e+01, 5.3791050276885002e+01, 1.1039779227737998e+01);
  mol.emplace_back(AtomicNumber(6), 4.0876662868059000e+01, 5.0346079798938000e+01, 7.0505676649589999e+00);
  mol.emplace_back(AtomicNumber(6), 3.8191362237690001e+01, 4.9932229807346999e+01, 5.9998800150750000e+00);
  mol.emplace_back(AtomicNumber(6), 3.8300966345052004e+01, 5.0245924321520995e+01, 3.1293862377839998e+00);
  mol.emplace_back(AtomicNumber(6), 3.6286518440778003e+01, 4.8862644897572999e+01, 1.6837458561990000e+00);
  mol.emplace_back(AtomicNumber(7), 3.3795859587275999e+01, 5.0160886652016003e+01, 2.0314554381750001e+00);
  mol.emplace_back(AtomicNumber(7), 4.1761054630910998e+01, 5.5110079017206999e+01, 1.0591914168345001e+01);
  mol.emplace_back(AtomicNumber(6), 4.1398227241023001e+01, 5.7755695401806996e+01, 1.1113478541309000e+01);
  mol.emplace_back(AtomicNumber(6), 4.0564858079874000e+01, 5.8492688537516997e+01, 1.3721300406129000e+01);
  mol.emplace_back(AtomicNumber(8), 3.9808967684273995e+01, 6.0682880958768003e+01, 1.4235305875137001e+01);
  mol.emplace_back(AtomicNumber(6), 3.9727709466747001e+01, 5.9002914554546997e+01, 9.0404491313760005e+00);
  mol.emplace_back(AtomicNumber(6), 4.1311299845529000e+01, 5.9227791947237996e+01, 6.5649080857860005e+00);
  mol.emplace_back(AtomicNumber(6), 3.9979043023284000e+01, 5.8063720738013998e+01, 4.3671567605790003e+00);
  mol.emplace_back(AtomicNumber(8), 3.7684915672637999e+01, 5.8190332379276995e+01, 4.1007053961299995e+00);
  mol.emplace_back(AtomicNumber(8), 4.1487044362505998e+01, 5.6979018020327999e+01, 2.9007293931149998e+00);
  mol.emplace_back(AtomicNumber(7), 4.0957921085586001e+01, 5.6756030353625995e+01, 1.5478745575899001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0476040958391003e+01, 5.7169880345217003e+01, 1.8179164014179999e+01);
  mol.emplace_back(AtomicNumber(6), 4.2526393656456001e+01, 5.9012363184492003e+01, 1.9154262624504000e+01);
  mol.emplace_back(AtomicNumber(8), 4.4557849094631003e+01, 5.9188107701469001e+01, 1.8054442098906001e+01);
  mol.emplace_back(AtomicNumber(6), 4.0888001223993001e+01, 5.4656544779846996e+01, 1.9564333164116999e+01);
  mol.emplace_back(AtomicNumber(8), 3.8822530718015997e+01, 5.3001144813483002e+01, 1.9008753723350999e+01);
  mol.emplace_back(AtomicNumber(7), 4.2029395721348997e+01, 6.0231236447397002e+01, 2.1242409842348998e+01);
  mol.emplace_back(AtomicNumber(6), 4.3864319656668002e+01, 6.1911202851618000e+01, 2.2470731735198999e+01);
  mol.emplace_back(AtomicNumber(6), 4.4425568275400998e+01, 6.0894530269535991e+01, 2.5114458393809997e+01);
  mol.emplace_back(AtomicNumber(8), 4.2601982696016002e+01, 6.0361627540637997e+01, 2.6520414529626002e+01);
  mol.emplace_back(AtomicNumber(6), 4.2894890224311006e+01, 6.4755240465063011e+01, 2.2648365978165000e+01);
  mol.emplace_back(AtomicNumber(8), 4.2509386122555000e+01, 6.5554594558410002e+01, 2.0010308497521002e+01);
  mol.emplace_back(AtomicNumber(6), 4.4837528541003003e+01, 6.6387963719558996e+01, 2.4041094032057998e+01);
  mol.emplace_back(AtomicNumber(7), 4.6846307267309996e+01, 6.0510915893769003e+01, 2.5734288518202000e+01);
  mol.emplace_back(AtomicNumber(6), 4.7524718897361005e+01, 5.9732348786301003e+01, 2.8308095315220001e+01);
  mol.emplace_back(AtomicNumber(6), 4.8562178465321999e+01, 6.2126631614363994e+01, 2.9610116521641000e+01);
  mol.emplace_back(AtomicNumber(8), 4.9431452420262005e+01, 6.3740457608969997e+01, 2.8145578880165999e+01);
  mol.emplace_back(AtomicNumber(6), 4.9718690770589994e+01, 5.7814276907466002e+01, 2.8283528877363000e+01);
  mol.emplace_back(AtomicNumber(6), 4.9680896250810001e+01, 5.5709122155720003e+01, 2.6380574806440002e+01);
  mol.emplace_back(AtomicNumber(6), 5.1765264016677001e+01, 5.3747586579138002e+01, 2.6888911097481000e+01);
  mol.emplace_back(AtomicNumber(6), 4.7133545617637999e+01, 5.4437336565122997e+01, 2.6365456998528000e+01);
  mol.emplace_back(AtomicNumber(7), 4.8416669564168998e+01, 6.2257022707605003e+01, 3.2030855513550001e+01);
  mol.emplace_back(AtomicNumber(6), 4.9471136666030993e+01, 6.4490678826603002e+01, 3.3353663705849996e+01);
  mol.emplace_back(AtomicNumber(6), 5.1920221547775000e+01, 6.3591169255839006e+01, 3.4589544502655997e+01);
  mol.emplace_back(AtomicNumber(8), 5.1980692779423002e+01, 6.1580500803543003e+01, 3.5825425299461997e+01);
  mol.emplace_back(AtomicNumber(6), 4.7647551086645997e+01, 6.5318378809784988e+01, 3.5489054073420000e+01);
  mol.emplace_back(AtomicNumber(6), 4.5311849764242005e+01, 6.6369066459669000e+01, 3.4253173276614000e+01);
  mol.emplace_back(AtomicNumber(7), 4.5075634015617005e+01, 6.8846497231248009e+01, 3.3601217810408997e+01);
  mol.emplace_back(AtomicNumber(6), 4.3131105972936005e+01, 6.5222002784346003e+01, 3.3603107536397999e+01);
  mol.emplace_back(AtomicNumber(6), 4.2847647074586000e+01, 6.9214993799102999e+01, 3.2503287010800001e+01);
  mol.emplace_back(AtomicNumber(7), 4.1659009427505005e+01, 6.7000234939994996e+01, 3.2452264409096998e+01);
  mol.emplace_back(AtomicNumber(7), 5.3904433836224996e+01, 6.5095391143083006e+01, 3.4372226013921001e+01);
  mol.emplace_back(AtomicNumber(6), 5.6315724198188995e+01, 6.4524693894405004e+01, 3.5581650646881002e+01);
  mol.emplace_back(AtomicNumber(6), 5.6790045421427997e+01, 6.6219778106538001e+01, 3.7802078683956005e+01);
  mol.emplace_back(AtomicNumber(8), 5.6890200898845002e+01, 6.8606502030645004e+01, 3.7393897870331998e+01);
  mol.emplace_back(AtomicNumber(6), 5.8439776209824998e+01, 6.4825160326656004e+01, 3.3548305482716998e+01);
  mol.emplace_back(AtomicNumber(6), 6.1123187114204995e+01, 6.4596503481987000e+01, 3.4691589706062004e+01);
  mol.emplace_back(AtomicNumber(6), 6.1520029571895002e+01, 6.1950887097387003e+01, 3.5659129412430005e+01);
  mol.emplace_back(AtomicNumber(6), 6.3043148719028991e+01, 6.5178539086599002e+01, 3.2588324680305000e+01);
  mol.emplace_back(AtomicNumber(7), 5.6926105692636000e+01, 6.5257907578137008e+01, 4.0045183432898995e+01);
  mol.emplace_back(AtomicNumber(6), 5.7596958418730999e+01, 6.6837718504940995e+01, 4.2280729277885996e+01);
  mol.emplace_back(AtomicNumber(6), 6.0284148775089001e+01, 6.5970334275989998e+01, 4.2949692277992000e+01);
  mol.emplace_back(AtomicNumber(8), 6.0830279585909999e+01, 6.3676206925343998e+01, 4.2773947761015002e+01);
  mol.emplace_back(AtomicNumber(6), 5.5694004347808004e+01, 6.6482450019008994e+01, 4.4404781289521999e+01);
  mol.emplace_back(AtomicNumber(6), 5.6353518717969003e+01, 6.7948877386473001e+01, 4.6799064117585004e+01);
  mol.emplace_back(AtomicNumber(6), 5.3004924265461000e+01, 6.6998345214006008e+01, 4.3597868292219005e+01);
  mol.emplace_back(AtomicNumber(7), 6.1913092577606996e+01, 6.7710771911859013e+01, 4.3633773086010002e+01);
  mol.emplace_back(AtomicNumber(6), 6.4524693894405004e+01, 6.7032360281807996e+01, 4.4372655947709006e+01);
  mol.emplace_back(AtomicNumber(6), 6.4702328137370998e+01, 6.6807482889116997e+01, 4.7203465479230999e+01);
  mol.emplace_back(AtomicNumber(8), 6.3696993911223004e+01, 6.8402411623833004e+01, 4.8618870244992003e+01);
  mol.emplace_back(AtomicNumber(6), 6.6355838377745997e+01, 6.9095941061795997e+01, 4.3287953230023000e+01);
  mol.emplace_back(AtomicNumber(6), 6.7890295880814008e+01, 6.7990451358230999e+01, 4.1076973822892995e+01);
  mol.emplace_back(AtomicNumber(6), 6.6146078792967003e+01, 6.6299146598076007e+01, 3.9533067689880006e+01);
  mol.emplace_back(AtomicNumber(6), 6.9037359556137005e+01, 7.0084267754043012e+01, 3.9527398511913006e+01);
  mol.emplace_back(AtomicNumber(7), 6.6008128795770006e+01, 6.4976338405776005e+01, 4.8095416146039000e+01);
  mol.emplace_back(AtomicNumber(6), 6.6444655499229000e+01, 6.4579495948086006e+01, 5.0826070200144002e+01);
  mol.emplace_back(AtomicNumber(6), 6.9298141742618995e+01, 6.4810042518743998e+01, 5.1190787316020995e+01);
  mol.emplace_back(AtomicNumber(8), 7.0496228019645002e+01, 6.2801263792436991e+01, 5.0635207875255006e+01);
  mol.emplace_back(AtomicNumber(6), 6.5605617160112999e+01, 6.1907423399639995e+01, 5.1563063335854004e+01);
  mol.emplace_back(AtomicNumber(6), 6.7561483558728000e+01, 6.0573276851406007e+01, 5.3214683850240000e+01);
  mol.emplace_back(AtomicNumber(6), 6.7296921920268005e+01, 5.7782151565653002e+01, 5.2995475635516001e+01);
  mol.emplace_back(AtomicNumber(7), 6.6215998654559996e+01, 5.7167990619228000e+01, 5.0512375685969999e+01);
  mol.emplace_back(AtomicNumber(6), 6.4889411010282004e+01, 5.4996695457867006e+01, 5.0361197606849998e+01);
  mol.emplace_back(AtomicNumber(7), 6.4458553484790002e+01, 5.3738137949193003e+01, 5.2473911262552001e+01);
  mol.emplace_back(AtomicNumber(7), 6.4277139789846004e+01, 5.4153877666772999e+01, 4.8106754501973001e+01);
  mol.emplace_back(AtomicNumber(7), 7.0292137612833002e+01, 6.6890630832632993e+01, 5.1992031135357003e+01);
  mol.emplace_back(AtomicNumber(6), 7.3071924542651999e+01, 6.7089052061478000e+01, 5.2307615375520001e+01);
  mol.emplace_back(AtomicNumber(6), 7.3842932746163996e+01, 6.6010018521758994e+01, 5.4860635186658996e+01);
  mol.emplace_back(AtomicNumber(8), 7.2370836200732995e+01, 6.6038364411594003e+01, 5.6684220766044000e+01);
  mol.emplace_back(AtomicNumber(6), 7.3850491650119992e+01, 6.9808367759649002e+01, 5.1789830454533998e+01);
  mol.emplace_back(AtomicNumber(6), 7.4647956017478009e+01, 7.0562368429260005e+01, 4.9136655165977999e+01);
  mol.emplace_back(AtomicNumber(6), 7.3102160158475996e+01, 6.9252788318882992e+01, 4.7097640823846994e+01);
  mol.emplace_back(AtomicNumber(6), 7.4336151229293009e+01, 7.3423413576605995e+01, 4.8872093527517997e+01);
  mol.emplace_back(AtomicNumber(7), 7.6144619000765999e+01, 6.5029250733468004e+01, 5.4887091350505003e+01);
  mol.emplace_back(AtomicNumber(6), 7.7238770348396997e+01, 6.3876517880178000e+01, 5.7169880345217003e+01);
  mol.emplace_back(AtomicNumber(6), 7.8924405930584996e+01, 6.5817266470880995e+01, 5.8475681003615996e+01);
  mol.emplace_back(AtomicNumber(8), 8.1154282597605004e+01, 6.6129071259065995e+01, 5.7793489921586996e+01);
  mol.emplace_back(AtomicNumber(6), 7.8708977167839009e+01, 6.1470896696181008e+01, 5.6546270768846995e+01);
  mol.emplace_back(AtomicNumber(6), 7.8627718950311987e+01, 5.9420543998115996e+01, 5.8560718673121002e+01);
  mol.emplace_back(AtomicNumber(6), 7.9171960035143996e+01, 5.6842957749119996e+01, 5.7553494720983998e+01);
  mol.emplace_back(AtomicNumber(7), 8.1845922309578995e+01, 5.6191002282915001e+01, 5.7755695401806996e+01);
  mol.emplace_back(AtomicNumber(6), 8.3476755838085992e+01, 5.6512255701045000e+01, 5.5848961878905996e+01);
  mol.emplace_back(AtomicNumber(7), 8.2683070922705994e+01, 5.7281374178568001e+01, 5.3585070144084000e+01);
  mol.emplace_back(AtomicNumber(7), 8.5939068801752995e+01, 5.6173994749014000e+01, 5.6243914610607000e+01);
  mol.emplace_back(AtomicNumber(7), 7.7790570337185002e+01, 6.7143854115159002e+01, 6.0278479597122001e+01);
  mol.emplace_back(AtomicNumber(6), 7.9075584009704997e+01, 6.9069484897949991e+01, 6.1767583676454002e+01);
  mol.emplace_back(AtomicNumber(6), 7.7953086772238990e+01, 7.1698093748649001e+01, 6.1582390529531999e+01);
  mol.emplace_back(AtomicNumber(8), 7.7671517599877987e+01, 7.2797914274247006e+01, 5.9526368653500001e+01);
  mol.emplace_back(AtomicNumber(7), 7.7376720345593995e+01, 7.2701538248808006e+01, 6.3791480210672994e+01);
  mol.emplace_back(AtomicNumber(6), 7.6293907353896998e+01, 7.5235660800057005e+01, 6.4144858970615999e+01);
  mol.emplace_back(AtomicNumber(6), 7.5647621065658996e+01, 7.5573921752087998e+01, 6.6956771242248010e+01);
  mol.emplace_back(AtomicNumber(8), 7.3572701929736994e+01, 7.6581145704224994e+01, 6.7438651369442994e+01);
  mol.emplace_back(AtomicNumber(8), 7.7217983362517998e+01, 7.4785906014675007e+01, 6.8504456827238997e+01);
  mol.emplace_back(AtomicNumber(1), 5.0878207740180514e+01, 4.4208344272965448e+01, 5.1067955126735995e+00);
  mol.emplace_back(AtomicNumber(1), 5.2471341235206957e+01, 4.6439052419420612e+01, 3.0019620143457302e+00);
  mol.emplace_back(AtomicNumber(1), 4.8762149269477980e+01, 4.8549649582014929e+01, 3.4728628335446401e+00);
  mol.emplace_back(AtomicNumber(1), 4.5732616152952737e+01, 4.8251148464792493e+01, 6.5078383609181998e+00);
  mol.emplace_back(AtomicNumber(1), 4.7169866151146579e+01, 4.4969147956136851e+01, 6.2884600708551899e+00);
  mol.emplace_back(AtomicNumber(1), 4.9755464838335939e+01, 4.5955698305954186e+01, 1.0137453965250391e+01);
  mol.emplace_back(AtomicNumber(1), 4.8062364838491384e+01, 4.8981149614343188e+01, 1.0481327403468720e+01);
  mol.emplace_back(AtomicNumber(1), 4.6413730091168013e+01, 4.3302693092737201e+01, 1.5115313473694520e+01);
  mol.emplace_back(AtomicNumber(1), 4.8085287214737960e+01, 4.6348818003445864e+01, 1.4565989025952110e+01);
  mol.emplace_back(AtomicNumber(1), 4.4716416002367993e+01, 4.6348345571948606e+01, 1.5580620703965991e+01);
  mol.emplace_back(AtomicNumber(1), 4.8175956267690182e+01, 5.2547572781603222e+01, 4.8715435243029903e+00);
  mol.emplace_back(AtomicNumber(1), 5.2823491673257109e+01, 5.4694566066745679e+01, 7.8903807917304905e+00);
  mol.emplace_back(AtomicNumber(1), 5.2036212928979822e+01, 5.6665531376012787e+01, 3.9921406380619504e+00);
  mol.emplace_back(AtomicNumber(1), 4.8481959597088952e+01, 5.6865313207569876e+01, 4.7939324779347601e+00);
  mol.emplace_back(AtomicNumber(1), 4.9238738163903783e+01, 6.0013313246345518e+01, 7.8719748605976312e+00);
  mol.emplace_back(AtomicNumber(1), 5.2740532702340005e+01, 5.9801474962978617e+01, 7.3577237272110594e+00);
  mol.emplace_back(AtomicNumber(1), 4.6748702919978150e+01, 6.0802973045368951e+01, 4.2204762293128200e+00);
  mol.emplace_back(AtomicNumber(1), 4.7964949463758437e+01, 6.3174805928682630e+01, 1.8862110986604601e+00);
  mol.emplace_back(AtomicNumber(1), 5.2771070674322246e+01, 5.6107722058579775e+01, 1.1539082628551581e+01);
  mol.emplace_back(AtomicNumber(1), 4.7443725241472457e+01, 5.6923119925573381e+01, 1.3892547375252180e+01);
  mol.emplace_back(AtomicNumber(1), 4.8957584731260361e+01, 5.5738582983888513e+01, 1.8139649843750011e+01);
  mol.emplace_back(AtomicNumber(1), 5.3357339265149605e+01, 5.2982512115231458e+01, 1.5525818650284988e+01);
  mol.emplace_back(AtomicNumber(1), 5.3614814431150862e+01, 5.6200129659441870e+01, 1.7039356886654762e+01);
  mol.emplace_back(AtomicNumber(1), 4.6073730591227132e+01, 5.2953674896639320e+01, 1.5792855829790581e+01);
  mol.emplace_back(AtomicNumber(1), 4.8637729710362223e+01, 5.1843479775361708e+01, 1.3655537941711801e+01);
  mol.emplace_back(AtomicNumber(1), 4.8637502943243540e+01, 5.0895649911058982e+01, 1.7043816639988798e+01);
  mol.emplace_back(AtomicNumber(1), 5.0901205705466637e+01, 5.2688754210241406e+01, 2.0517586542008161e+01);
  mol.emplace_back(AtomicNumber(1), 5.3932288397302855e+01, 5.1223555164670259e+01, 1.9497569144925631e+01);
  mol.emplace_back(AtomicNumber(1), 5.3931778171285835e+01, 5.4412581154667102e+01, 2.0983800840754352e+01);
  mol.emplace_back(AtomicNumber(1), 4.7809670679242316e+01, 5.9915916768872464e+01, 1.7383740550890121e+01);
  mol.emplace_back(AtomicNumber(1), 5.2355671107420271e+01, 6.3706272465828995e+01, 1.6258352032660952e+01);
  mol.emplace_back(AtomicNumber(1), 4.6814616562474470e+01, 6.4685188322650774e+01, 1.8067443413710318e+01);
  mol.emplace_back(AtomicNumber(1), 4.9393601208702329e+01, 6.7092623643597207e+01, 1.7477319781865400e+01);
  mol.emplace_back(AtomicNumber(1), 4.4421524261784540e+01, 6.3129433607686742e+01, 1.4623587874096831e+01);
  mol.emplace_back(AtomicNumber(1), 5.1067633873317874e+01, 6.7781201999469033e+01, 1.3285265031427141e+01);
  mol.emplace_back(AtomicNumber(1), 4.3333042092120536e+01, 6.3307880432828007e+01, 9.9814003930787703e+00);
  mol.emplace_back(AtomicNumber(1), 5.0022369737022302e+01, 6.7831827758714340e+01, 8.5856298803434790e+00);
  mol.emplace_back(AtomicNumber(1), 4.5945304813014694e+01, 6.5954573961241735e+01, 7.1853807170142598e+00);
  mol.emplace_back(AtomicNumber(1), 5.4663461176966734e+01, 6.5034447479937739e+01, 1.9625730361499610e+01);
  mol.emplace_back(AtomicNumber(1), 5.2609650280341867e+01, 6.3036931520525187e+01, 2.4737722620642963e+01);
  mol.emplace_back(AtomicNumber(1), 5.8188215886169317e+01, 6.4002864959802537e+01, 2.3166188693670779e+01);
  mol.emplace_back(AtomicNumber(1), 5.7225400494773822e+01, 6.1114040840418234e+01, 2.7841068434298542e+01);
  mol.emplace_back(AtomicNumber(1), 5.8934827727163331e+01, 6.4156461888188460e+01, 2.7398022176177488e+01);
  mol.emplace_back(AtomicNumber(1), 5.5464156975765931e+01, 6.4155970559431324e+01, 2.7975371260336768e+01);
  mol.emplace_back(AtomicNumber(1), 5.4577025110229876e+01, 5.9486193078973855e+01, 2.3129830365642420e+01);
  mol.emplace_back(AtomicNumber(1), 5.7623452377096783e+01, 6.0000614287699442e+01, 2.1447728571053851e+01);
  mol.emplace_back(AtomicNumber(1), 5.7622979945599532e+01, 5.8951891952844001e+01, 2.4806149598704650e+01);
  mol.emplace_back(AtomicNumber(1), 5.1105995310894571e+01, 6.6446923170415801e+01, 2.6756271230313093e+01);
  mol.emplace_back(AtomicNumber(1), 5.3197884186197790e+01, 7.1827369903556502e+01, 2.5255696514227857e+01);
  mol.emplace_back(AtomicNumber(1), 4.8545265417720451e+01, 7.0696104337501524e+01, 2.5200157467411149e+01);
  mol.emplace_back(AtomicNumber(1), 4.8902858266618921e+01, 7.0123838616252655e+01, 2.8780262250831541e+01);
  mol.emplace_back(AtomicNumber(1), 4.9732466873049809e+01, 7.4505036040929923e+01, 2.9535529036855170e+01);
  mol.emplace_back(AtomicNumber(1), 4.9758261632799659e+01, 7.5246847977911884e+01, 2.6109455818798171e+01);
  mol.emplace_back(AtomicNumber(1), 4.5216758752475521e+01, 7.3414267302819241e+01, 2.5853586919887569e+01);
  mol.emplace_back(AtomicNumber(1), 4.5323641654413358e+01, 7.3245628155560880e+01, 2.9441987600399671e+01);
  mol.emplace_back(AtomicNumber(1), 4.6005209126865985e+01, 7.7899550834970626e+01, 2.9582394241382371e+01);
  mol.emplace_back(AtomicNumber(1), 4.5642098278079644e+01, 7.8045513270360985e+01, 2.6094867134163088e+01);
  mol.emplace_back(AtomicNumber(1), 4.1584025100261485e+01, 7.4953146764701501e+01, 2.8244222576791799e+01);
  mol.emplace_back(AtomicNumber(1), 4.1761545959668140e+01, 7.7946132580599482e+01, 2.9963192925425762e+01);
  mol.emplace_back(AtomicNumber(1), 5.5944241863271380e+01, 7.3773712083186922e+01, 2.7531947057017920e+01);
  mol.emplace_back(AtomicNumber(1), 5.7001203403438858e+01, 7.1264288250614172e+01, 3.2722967657021250e+01);
  mol.emplace_back(AtomicNumber(1), 6.1267675563323948e+01, 7.3055880769005384e+01, 3.2482896867378692e+01);
  mol.emplace_back(AtomicNumber(1), 5.8011450917158264e+01, 7.6721798009586266e+01, 2.9679960794194439e+01);
  mol.emplace_back(AtomicNumber(1), 5.9290663130892028e+01, 7.0897964867646508e+01, 2.7611636801974051e+01);
  mol.emplace_back(AtomicNumber(1), 6.2160212045188537e+01, 7.2909313621298551e+01, 2.7303422493168149e+01);
  mol.emplace_back(AtomicNumber(1), 6.2159645127391833e+01, 7.0171062868717769e+01, 2.9512663352388270e+01);
  mol.emplace_back(AtomicNumber(1), 5.8419518347222926e+01, 7.4521155403616092e+01, 3.5978852152508914e+01);
  mol.emplace_back(AtomicNumber(1), 5.3807679865588199e+01, 7.7927273115229269e+01, 3.6578896845796081e+01);
  mol.emplace_back(AtomicNumber(1), 5.8650367274039155e+01, 7.6872106814751334e+01, 3.9718185247762435e+01);
  mol.emplace_back(AtomicNumber(1), 5.5867972522355345e+01, 7.9015604106814138e+01, 4.0638103859207639e+01);
  mol.emplace_back(AtomicNumber(1), 5.6073215662020630e+01, 7.3280833750735951e+01, 4.0003628358400888e+01);
  mol.emplace_back(AtomicNumber(1), 5.5446336859689659e+01, 7.3099476747571629e+01, 4.4397071207486881e+01);
  mol.emplace_back(AtomicNumber(1), 5.7194616858413006e+01, 7.6145393788421487e+01, 4.4195116191042445e+01);
  mol.emplace_back(AtomicNumber(1), 5.3688154696783954e+01, 7.6144902459664351e+01, 4.4484282061879227e+01);
  mol.emplace_back(AtomicNumber(1), 5.1647137245104609e+01, 7.3030615132532461e+01, 3.9948882996499556e+01);
  mol.emplace_back(AtomicNumber(1), 5.1239126506819623e+01, 7.6077061296659238e+01, 4.1659878701459938e+01);
  mol.emplace_back(AtomicNumber(1), 5.2054165325875317e+01, 7.6076569967902103e+01, 3.8237244784702916e+01);
  mol.emplace_back(AtomicNumber(1), 5.9321503459032513e+01, 7.8742406420584402e+01, 3.3172193319126329e+01);
  mol.emplace_back(AtomicNumber(1), 5.8771498709934058e+01, 8.4524892357884852e+01, 3.4602432433900979e+01);
  mol.emplace_back(AtomicNumber(1), 6.3183744332610601e+01, 8.4961834801061428e+01, 3.2822140476923970e+01);
  mol.emplace_back(AtomicNumber(1), 6.0839841599414335e+01, 8.0830289076790947e+01, 2.9967463706160899e+01);
  mol.emplace_back(AtomicNumber(1), 6.3536423893937673e+01, 7.9711098859805702e+01, 3.5165362908764187e+01);
  mol.emplace_back(AtomicNumber(1), 6.5654353193369303e+01, 8.2513014480955889e+01, 3.5361346391083380e+01);
  mol.emplace_back(AtomicNumber(1), 6.2574610057316335e+01, 8.2512712124797645e+01, 3.7062534418160844e+01);
  mol.emplace_back(AtomicNumber(1), 5.5547966323378077e+01, 8.0960850245370963e+01, 3.1176812750081339e+01);
  mol.emplace_back(AtomicNumber(1), 5.2837740207214168e+01, 8.2115378338350510e+01, 2.7905489193263548e+01);
  mol.emplace_back(AtomicNumber(1), 5.4672323991855144e+01, 8.5221105206752227e+01, 2.7610748630759218e+01);
  mol.emplace_back(AtomicNumber(1), 5.8397692012049966e+01, 8.0000850545699066e+01, 2.7366104704223279e+01);
  mol.emplace_back(AtomicNumber(1), 5.8795214771096006e+01, 8.0758233824830370e+01, 2.1643768745152713e+01);
  mol.emplace_back(AtomicNumber(1), 6.2775109984749243e+01, 8.0673592997783061e+01, 2.4062977059010620e+01);
  mol.emplace_back(AtomicNumber(1), 6.1964474227247912e+01, 7.7429519289726656e+01, 2.5249460418464160e+01);
  mol.emplace_back(AtomicNumber(1), 6.1863789626553995e+01, 7.5802937644694921e+01, 2.0759698235718840e+01);
  mol.emplace_back(AtomicNumber(1), 6.2674368692275650e+01, 7.9001525648196093e+01, 1.9559400979285712e+01);
  mol.emplace_back(AtomicNumber(1), 6.6439118602081237e+01, 7.8560784855781620e+01, 2.3069340236734529e+01);
  mol.emplace_back(AtomicNumber(1), 6.5809008368309065e+01, 7.5082290638789758e+01, 2.2423564174513558e+01);
  mol.emplace_back(AtomicNumber(1), 6.7475406539929040e+01, 7.9247360102105091e+01, 1.8737445763110273e+01);
  mol.emplace_back(AtomicNumber(1), 6.9186685703787788e+01, 7.6296874223699731e+01, 1.9886550342501391e+01);
  mol.emplace_back(AtomicNumber(1), 6.4213191462378234e+01, 7.5747096241719959e+01, 1.7227970437616850e+01);
  mol.emplace_back(AtomicNumber(1), 6.7017677110873464e+01, 7.6588817991740342e+01, 1.5392007153003899e+01);
  mol.emplace_back(AtomicNumber(1), 5.7388748409262980e+01, 7.7503067425218546e+01, 1.9020129873804780e+01);
  mol.emplace_back(AtomicNumber(1), 5.5532470570268281e+01, 7.2804131472750811e+01, 2.2050683442364083e+01);
  mol.emplace_back(AtomicNumber(1), 5.3130118612232252e+01, 7.4807316610130357e+01, 1.7025637475974619e+01);
  mol.emplace_back(AtomicNumber(1), 5.0566006109537817e+01, 7.3648857887093698e+01, 2.1940285650086700e+01);
  mol.emplace_back(AtomicNumber(1), 4.9979548546111559e+01, 7.1136807335396213e+01, 1.8362278462514102e+01);
  mol.emplace_back(AtomicNumber(1), 5.2843409385181168e+01, 7.0485683348626381e+01, 1.6426027419664919e+01);
  mol.emplace_back(AtomicNumber(1), 5.2843069234503147e+01, 6.9755927863454247e+01, 1.9867860952470178e+01);
  mol.emplace_back(AtomicNumber(1), 5.7442133168452230e+01, 6.9092766322134480e+01, 2.0776006571003908e+01);
  mol.emplace_back(AtomicNumber(1), 6.0695390944815180e+01, 6.9819951779961571e+01, 1.5750563762156760e+01);
  mol.emplace_back(AtomicNumber(1), 6.3867711757069259e+01, 6.8900297730154833e+01, 1.7699589352691579e+01);
  mol.emplace_back(AtomicNumber(1), 6.1231676283233497e+01, 6.3847794045145200e+01, 1.9276395615173069e+01);
  mol.emplace_back(AtomicNumber(1), 6.3380238040946821e+01, 6.4474276005018481e+01, 1.6545495896689498e+01);
  mol.emplace_back(AtomicNumber(1), 6.0593175666070167e+01, 6.9021958289326648e+01, 2.2428364078525618e+01);
  mol.emplace_back(AtomicNumber(1), 6.3631647186523388e+01, 6.7283108023288406e+01, 2.2773163482478562e+01);
  mol.emplace_back(AtomicNumber(1), 6.3631136960506353e+01, 7.0622367229410742e+01, 2.1664933676229509e+01);
  mol.emplace_back(AtomicNumber(1), 6.5455572916586405e+01, 6.2162555305414891e+01, 2.0396681873231937e+01);
  mol.emplace_back(AtomicNumber(1), 6.7041884500792563e+01, 6.5129462902664670e+01, 1.9369143366713189e+01);
  mol.emplace_back(AtomicNumber(1), 6.4993213658857769e+01, 6.5129066060206981e+01, 2.2229527109963040e+01);
  mol.emplace_back(AtomicNumber(1), 6.0447175436160030e+01, 6.7227606770991471e+01, 1.2614676866970600e+01);
  mol.emplace_back(AtomicNumber(1), 5.5858240433511995e+01, 6.3468941778870480e+01, 1.3368016132485449e+01);
  mol.emplace_back(AtomicNumber(1), 5.5452780825312146e+01, 6.3733597903629935e+01, 8.5649184835040391e+00);
  mol.emplace_back(AtomicNumber(1), 6.0192213605724149e+01, 6.3323470672237264e+01, 8.1378215127301505e+00);
  mol.emplace_back(AtomicNumber(1), 5.3456171934374311e+01, 6.7317784495186572e+01, 9.6659484337350001e+00);
  mol.emplace_back(AtomicNumber(1), 5.6190151906219945e+01, 6.8608202784035100e+01, 1.1464532938285529e+01);
  mol.emplace_back(AtomicNumber(1), 5.6189887344581493e+01, 6.8713479418882287e+01, 7.9477528727565305e+00);
  mol.emplace_back(AtomicNumber(1), 5.6359528046614017e+01, 5.9356898026806476e+01, 1.2988842612792601e+01);
  mol.emplace_back(AtomicNumber(1), 6.1899203091587857e+01, 5.7585298809378870e+01, 1.1957373476216729e+01);
  mol.emplace_back(AtomicNumber(1), 5.8238236933098150e+01, 5.7386424046296511e+01, 1.6206384567963450e+01);
  mol.emplace_back(AtomicNumber(1), 5.9297692911571112e+01, 5.4048015216869217e+01, 1.4639631647743439e+01);
  mol.emplace_back(AtomicNumber(1), 6.2577180084661386e+01, 5.3747359812019319e+01, 1.6651641805491629e+01);
  mol.emplace_back(AtomicNumber(1), 6.5728317068578775e+01, 5.6291611294569357e+01, 1.7017832907640049e+01);
  mol.emplace_back(AtomicNumber(1), 6.3781105614993393e+01, 5.8670379472262667e+01, 1.5307687579374718e+01);
  mol.emplace_back(AtomicNumber(1), 6.3870338476193965e+01, 5.8669793657206078e+01, 1.8824902281881190e+01);
  mol.emplace_back(AtomicNumber(1), 6.0043511067649739e+01, 5.3622165465248074e+01, 2.0370660346363412e+01);
  mol.emplace_back(AtomicNumber(1), 6.1622679487617475e+01, 5.6668592732114966e+01, 2.1145448001853410e+01);
  mol.emplace_back(AtomicNumber(1), 5.8474093633785238e+01, 5.6668101403357824e+01, 1.9575350266632867e+01);
  mol.emplace_back(AtomicNumber(1), 6.2178315620163154e+01, 5.3781771722279011e+01, 9.9895073175715794e+00);
  mol.emplace_back(AtomicNumber(1), 5.6914238213425079e+01, 5.2041126216551220e+01, 7.6759724810185510e+00);
  mol.emplace_back(AtomicNumber(1), 6.2288581131621299e+01, 5.1608945882866919e+01, 5.6552317892011805e+00);
  mol.emplace_back(AtomicNumber(1), 5.9492126818579322e+01, 4.9665287114140860e+01, 4.5833414137205999e+00);
  mol.emplace_back(AtomicNumber(1), 5.7128136298119990e+01, 5.3848271179831919e+01, 3.9030022631608206e+00);
  mol.emplace_back(AtomicNumber(1), 6.0392826916716388e+01, 5.5404120381095403e+01, 4.3606750004367303e+00);
  mol.emplace_back(AtomicNumber(1), 6.3264416735081014e+01, 5.2207289822763990e+01, 2.0636563690275600e+00);
  mol.emplace_back(AtomicNumber(1), 5.5509453707722258e+01, 4.8874814732942163e+01, 9.4878039647519703e+00);
  mol.emplace_back(AtomicNumber(1), 5.9207628570935370e+01, 4.4950666435964429e+01, 1.1638009784075731e+01);
  mol.emplace_back(AtomicNumber(1), 5.6633557212278909e+01, 4.4617243182465273e+01, 1.5579581354672040e+01);
  mol.emplace_back(AtomicNumber(1), 5.9754723142010761e+01, 4.6991929552022334e+01, 1.6935516443559212e+01);
  mol.emplace_back(AtomicNumber(1), 5.9349962732426846e+01, 4.9423647851927434e+01, 1.4425979227427099e+01);
  mol.emplace_back(AtomicNumber(1), 5.7230068117966646e+01, 4.9423062036870839e+01, 1.7233998663521760e+01);
  mol.emplace_back(AtomicNumber(1), 5.2416274619887503e+01, 4.5857583732605313e+01, 1.4599456073217299e+01);
  mol.emplace_back(AtomicNumber(1), 5.3552888110491331e+01, 4.8703208715881068e+01, 1.6327251442219890e+01);
  mol.emplace_back(AtomicNumber(1), 5.3454452283724322e+01, 4.8702641798084372e+01, 1.2810282404091989e+01);
  mol.emplace_back(AtomicNumber(1), 5.7035615313698550e+01, 4.1406863288792735e+01, 1.3123032055271491e+01);
  mol.emplace_back(AtomicNumber(1), 5.2382089476746486e+01, 4.0127745561358410e+01, 9.7742486301645908e+00);
  mol.emplace_back(AtomicNumber(1), 5.6897155090484517e+01, 3.6923526174410014e+01, 1.1703696659453369e+01);
  mol.emplace_back(AtomicNumber(1), 5.3777406455244417e+01, 3.5609391824399516e+01, 1.0467041074991879e+01);
  mol.emplace_back(AtomicNumber(1), 5.7708830197279802e+01, 3.8845339710703229e+01, 7.6000243935206395e+00);
  mol.emplace_back(AtomicNumber(1), 5.6823380187873958e+01, 3.5319054323449556e+01, 7.2550738114885798e+00);
  mol.emplace_back(AtomicNumber(1), 5.7443040236926947e+01, 3.8692101830255218e+01, 3.5261720036943296e+00);
  mol.emplace_back(AtomicNumber(1), 4.8620684381941444e+01, 3.9083596363396353e+01, 1.7459008337031992e+01);
  mol.emplace_back(AtomicNumber(1), 4.6296529285330230e+01, 3.4423040745765213e+01, 1.4592086141860200e+01);
  mol.emplace_back(AtomicNumber(1), 4.4622477723454800e+01, 3.6782665999189945e+01, 1.6428087220992929e+01);
  mol.emplace_back(AtomicNumber(1), 4.4470449267639751e+01, 3.7139597443992272e+01, 1.1357423269229010e+01);
  mol.emplace_back(AtomicNumber(1), 4.5192759232415220e+01, 4.0042405535695174e+01, 1.3037030625512100e+01);
  mol.emplace_back(AtomicNumber(1), 4.8897019013312907e+01, 3.9661814721510574e+01, 1.0141422389827289e+01);
  mol.emplace_back(AtomicNumber(1), 4.8998138250984297e+01, 3.6266468448034708e+01, 1.0261382195609009e+01);
  mol.emplace_back(AtomicNumber(1), 5.2012402381518420e+01, 3.4306142296085667e+01, 1.4185662773405971e+01);
  mol.emplace_back(AtomicNumber(1), 5.1605865629504848e+01, 3.0893259365431888e+01, 1.8749785673818440e+01);
  mol.emplace_back(AtomicNumber(1), 5.5552879610949482e+01, 2.8784533031566681e+01, 1.6284486943088822e+01);
  mol.emplace_back(AtomicNumber(1), 5.2242797774097305e+01, 2.8599037528486441e+01, 1.5001551969156720e+01);
  mol.emplace_back(AtomicNumber(1), 5.3033402436115232e+01, 3.1644822288517076e+01, 1.2060968254933710e+01);
  mol.emplace_back(AtomicNumber(1), 5.4852584953945858e+01, 3.6395498938563627e+01, 1.7365391311536932e+01);
  mol.emplace_back(AtomicNumber(1), 5.9838948229340495e+01, 3.4704987863324007e+01, 1.9648747224045628e+01);
  mol.emplace_back(AtomicNumber(1), 5.7057422751611611e+01, 3.9879019826686225e+01, 1.8664615723494208e+01);
  mol.emplace_back(AtomicNumber(1), 6.0099257984325234e+01, 3.9750990890931476e+01, 2.0462387645869470e+01);
  mol.emplace_back(AtomicNumber(1), 6.1292563254599067e+01, 3.5029699480013875e+01, 1.7214251026936708e+01);
  mol.emplace_back(AtomicNumber(1), 6.1344001596019645e+01, 3.5233411941628077e+01, 2.3482925666687070e+01);
  mol.emplace_back(AtomicNumber(1), 5.7396987614575025e+01, 3.5341012939441740e+01, 2.7645897534154621e+01);
  mol.emplace_back(AtomicNumber(1), 6.1226819687441754e+01, 3.4469093368117136e+01, 3.0481828223106806e+01);
  mol.emplace_back(AtomicNumber(1), 6.5114156122153773e+01, 3.3847430209515814e+01, 2.6953388548225682e+01);
  mol.emplace_back(AtomicNumber(1), 5.9048834896079697e+01, 3.1041697341867838e+01, 2.7520363036705351e+01);
  mol.emplace_back(AtomicNumber(1), 6.2082827765938980e+01, 3.1192214016891690e+01, 2.5746382764531599e+01);
  mol.emplace_back(AtomicNumber(1), 6.2082374231701628e+01, 3.0357049616053143e+01, 2.9164178982756781e+01);
  mol.emplace_back(AtomicNumber(1), 5.7798837846135875e+01, 3.8004978563394928e+01, 3.1054906526530949e+01);
  mol.emplace_back(AtomicNumber(1), 5.8308515842629063e+01, 4.3424561521767806e+01, 2.9678278938064228e+01);
  mol.emplace_back(AtomicNumber(1), 5.8254412987563995e+01, 4.0690241399244151e+01, 3.5052413986401660e+01);
  mol.emplace_back(AtomicNumber(1), 5.4608848095884646e+01, 3.9613664503310851e+01, 3.2303996508000054e+01);
  mol.emplace_back(AtomicNumber(1), 5.3734509675294120e+01, 4.2087901637968443e+01, 3.4774378601640095e+01);
  mol.emplace_back(AtomicNumber(1), 5.6419394565945538e+01, 4.5808885493868779e+01, 3.5126491245170456e+01);
  mol.emplace_back(AtomicNumber(1), 5.9403347491616096e+01, 4.4556431800139251e+01, 3.6505651066462441e+01);
  mol.emplace_back(AtomicNumber(1), 5.9402818368339176e+01, 4.5931245251656527e+01, 3.3267038666514239e+01);
  mol.emplace_back(AtomicNumber(1), 5.1531353733758579e+01, 4.3278958134315360e+01, 3.1030075527035489e+01);
  mol.emplace_back(AtomicNumber(1), 5.4334535471321402e+01, 4.5392900111910215e+01, 3.1249378228058944e+01);
  mol.emplace_back(AtomicNumber(1), 5.4334233115163158e+01, 4.2700002783065436e+01, 2.8985089650779251e+01);
  mol.emplace_back(AtomicNumber(1), 6.2574364392937767e+01, 3.8853163176297691e+01, 3.2452434484436004e+01);
  mol.emplace_back(AtomicNumber(1), 6.6596665057784151e+01, 4.2540585498633391e+01, 3.3758027272976221e+01);
  mol.emplace_back(AtomicNumber(1), 6.6240281633518649e+01, 3.7731988747023991e+01, 3.4605815043421288e+01);
  mol.emplace_back(AtomicNumber(1), 6.7391275938898772e+01, 3.7450873108900346e+01, 3.1098823758515312e+01);
  mol.emplace_back(AtomicNumber(1), 7.1391599090493088e+01, 3.9231713086414167e+01, 3.2240917454487239e+01);
  mol.emplace_back(AtomicNumber(1), 7.0387058549260473e+01, 3.9398897144660999e+01, 3.5640326638839454e+01);
  mol.emplace_back(AtomicNumber(1), 7.0666435639474230e+01, 3.6575363058196650e+01, 3.7839797614696437e+01);
  mol.emplace_back(AtomicNumber(1), 6.4942153262634989e+01, 3.9447425308058520e+01, 2.7919718829960722e+01);
  mol.emplace_back(AtomicNumber(1), 6.9401339678878301e+01, 4.1281933503659936e+01, 2.5153991461499878e+01);
  mol.emplace_back(AtomicNumber(1), 6.7024045487456391e+01, 3.7630963995652053e+01, 2.3571345945712380e+01);
  mol.emplace_back(AtomicNumber(1), 6.4108765204226103e+01, 3.9539757319881055e+01, 2.3216323124158951e+01);
  mol.emplace_back(AtomicNumber(1), 6.3461175005055694e+01, 4.0467367116101492e+01, 1.9149254850633149e+01);
  mol.emplace_back(AtomicNumber(1), 6.6009583884781534e+01, 4.1422113377523964e+01, 1.6880374239200190e+01);
  mol.emplace_back(AtomicNumber(1), 6.3016673657923107e+01, 4.3522524917037572e+01, 2.5938397822273892e+01);
  mol.emplace_back(AtomicNumber(1), 6.3466220573446314e+01, 4.7815642213367546e+01, 2.2408295188522441e+01);
  mol.emplace_back(AtomicNumber(1), 6.0301231898029556e+01, 4.6986600524733355e+01, 2.7256746261239851e+01);
  mol.emplace_back(AtomicNumber(1), 5.7729201443441220e+01, 5.0749385119510379e+01, 2.5209757275435269e+01);
  mol.emplace_back(AtomicNumber(1), 6.0683617951903706e+01, 5.1518371316214150e+01, 2.6957583739921262e+01);
  mol.emplace_back(AtomicNumber(1), 6.0683221109446023e+01, 5.1466063700838632e+01, 2.3439613147019191e+01);
  mol.emplace_back(AtomicNumber(1), 5.9238147645657719e+01, 4.4057525241783360e+01, 2.3938954342352549e+01);
  mol.emplace_back(AtomicNumber(1), 5.6811115866205348e+01, 4.6602438128429547e+01, 2.4027034470699842e+01);
  mol.emplace_back(AtomicNumber(1), 5.9270084014871820e+01, 4.6602230258570764e+01, 2.1510637549227660e+01);
  mol.emplace_back(AtomicNumber(1), 6.4879508846099640e+01, 4.7104802885345308e+01, 2.9114309113907069e+01);
  mol.emplace_back(AtomicNumber(1), 6.6368783000770648e+01, 5.2442126071417022e+01, 2.9759990689828591e+01);
  mol.emplace_back(AtomicNumber(1), 6.7947119941303228e+01, 4.7410957392823200e+01, 3.2437373368303682e+01);
  mol.emplace_back(AtomicNumber(1), 6.9070108507526371e+01, 5.0558995226118633e+01, 3.3481069032028380e+01);
  mol.emplace_back(AtomicNumber(1), 6.4339009418725865e+01, 5.2047022161636903e+01, 3.3651219960077938e+01);
  mol.emplace_back(AtomicNumber(1), 6.3357618020858496e+01, 4.8767723961145528e+01, 3.2921804625583832e+01);
  mol.emplace_back(AtomicNumber(1), 6.5758136944685191e+01, 4.7199742719032670e+01, 3.6735347260425385e+01);
  mol.emplace_back(AtomicNumber(1), 6.6753512314871159e+01, 5.0514794535235922e+01, 3.7444013403560284e+01);
  mol.emplace_back(AtomicNumber(1), 6.2198403407426227e+01, 5.1800677481710864e+01, 3.7927140749908020e+01);
  mol.emplace_back(AtomicNumber(1), 6.1141063922060944e+01, 4.8570795615831834e+01, 3.7170966895409670e+01);
  mol.emplace_back(AtomicNumber(1), 6.1370627835204665e+01, 4.8684198072431727e+01, 4.1705402200534948e+01);
  mol.emplace_back(AtomicNumber(1), 6.4321907398525411e+01, 4.7155768795268635e+01, 4.0757610130751999e+01);
  mol.emplace_back(AtomicNumber(1), 6.9838225430275202e+01, 4.6881834115903203e+01, 2.7670104924073712e+01);
  mol.emplace_back(AtomicNumber(1), 7.4731878748649493e+01, 4.9597445951135761e+01, 2.7563505481034220e+01);
  mol.emplace_back(AtomicNumber(1), 7.2542536704093550e+01, 4.4690318886459899e+01, 2.5349993841078959e+01);
  mol.emplace_back(AtomicNumber(1), 7.5429508892008613e+01, 4.6018815153986786e+01, 2.3841557864879491e+01);
  mol.emplace_back(AtomicNumber(1), 7.5428941974211924e+01, 4.5261866511832949e+01, 2.7277533247118850e+01);
  mol.emplace_back(AtomicNumber(1), 6.9610381167781469e+01, 4.8877535938366321e+01, 2.2950174115868190e+01);
  mol.emplace_back(AtomicNumber(1), 7.2502266643267944e+01, 5.1365133435766140e+01, 1.8886866397060501e+01);
  mol.emplace_back(AtomicNumber(1), 6.6717550829300492e+01, 5.1446221577954134e+01, 2.0044002311904872e+01);
  mol.emplace_back(AtomicNumber(1), 6.8317808591305479e+01, 5.2417748606158916e+01, 1.7002110387411570e+01);
  mol.emplace_back(AtomicNumber(1), 6.9738788048734008e+01, 4.8180756171702242e+01, 1.6176961534314721e+01);
  mol.emplace_back(AtomicNumber(1), 6.8818756053729473e+01, 4.7088834700738261e+01, 1.9403177331775080e+01);
  mol.emplace_back(AtomicNumber(1), 6.4529437106637388e+01, 4.6958386915717590e+01, 1.8445709864928450e+01);
  mol.emplace_back(AtomicNumber(1), 6.4795397142329250e+01, 4.9582139170624856e+01, 1.6058721379182987e+01);
  mol.emplace_back(AtomicNumber(1), 6.7438783650262224e+01, 4.4544998957905797e+01, 1.4955556038584460e+01);
  mol.emplace_back(AtomicNumber(1), 6.3911912447951963e+01, 4.4727886639121223e+01, 1.4707019276511179e+01);
  mol.emplace_back(AtomicNumber(1), 6.4046517630148429e+01, 4.7532334493096670e+01, 1.1401453884772710e+01);
  mol.emplace_back(AtomicNumber(1), 6.7029601281864061e+01, 4.8958926436712552e+01, 1.2406674727361370e+01);
  mol.emplace_back(AtomicNumber(1), 6.8713195959983935e+01, 5.3631273944515051e+01, 2.4151643002414499e+01);
  mol.emplace_back(AtomicNumber(1), 6.8686002803002239e+01, 5.9008848294152465e+01, 2.2290225108729718e+01);
  mol.emplace_back(AtomicNumber(1), 6.7912670236523766e+01, 5.6596613069193957e+01, 2.7508627838313661e+01);
  mol.emplace_back(AtomicNumber(1), 6.4130440361319927e+01, 5.9617755802328041e+01, 2.4017963785952638e+01);
  mol.emplace_back(AtomicNumber(1), 6.4500448709966136e+01, 5.6119249387112674e+01, 2.3944604623059661e+01);
  mol.emplace_back(AtomicNumber(1), 6.5850450059247834e+01, 6.1560696475178283e+01, 2.7853257166927591e+01);
  mol.emplace_back(AtomicNumber(1), 6.8879076107298360e+01, 6.0459137401670404e+01, 2.9263484083478730e+01);
  mol.emplace_back(AtomicNumber(1), 6.8878565881281332e+01, 6.2078027861926913e+01, 2.6139710331882061e+01);
  mol.emplace_back(AtomicNumber(1), 6.2561759920591136e+01, 5.5561969192956568e+01, 2.7938011377534242e+01);
  mol.emplace_back(AtomicNumber(1), 6.3809394813048719e+01, 5.8598532090160887e+01, 2.9201897913497220e+01);
  mol.emplace_back(AtomicNumber(1), 6.1016020753368814e+01, 5.8598097453183421e+01, 2.7062784785728887e+01);
  mol.emplace_back(AtomicNumber(1), 7.2750765610821446e+01, 5.5040499306292020e+01, 2.6359863409600560e+01);
  mol.emplace_back(AtomicNumber(1), 7.5911313430164071e+01, 5.9369672574492114e+01, 2.7959988890786310e+01);
  mol.emplace_back(AtomicNumber(1), 7.6432669933269267e+01, 5.5135023400261801e+01, 2.9930651843895181e+01);
  mol.emplace_back(AtomicNumber(1), 7.7523608746718978e+01, 5.3782376434595491e+01, 2.6869201255415728e+01);
  mol.emplace_back(AtomicNumber(1), 8.1336603361023720e+01, 5.6549899042745878e+01, 2.6859714830950949e+01);
  mol.emplace_back(AtomicNumber(1), 8.0178371405105736e+01, 5.8328962677830042e+01, 2.9657605335744570e+01);
  mol.emplace_back(AtomicNumber(1), 8.3313710279755071e+01, 5.7279881295036688e+01, 3.2442456731214087e+01);
  mol.emplace_back(AtomicNumber(1), 8.3848030303144824e+01, 5.3933289952077033e+01, 3.3476212436236644e+01);
  mol.emplace_back(AtomicNumber(1), 7.6339789900909921e+01, 5.4798992324897824e+01, 2.2892650856763030e+01);
  mol.emplace_back(AtomicNumber(1), 8.0968579224665973e+01, 5.6827291920671186e+01, 2.0399554256735222e+01);
  mol.emplace_back(AtomicNumber(1), 7.7902064170535994e+01, 5.2684445634986488e+01, 1.9817461960343550e+01);
  mol.emplace_back(AtomicNumber(1), 7.6439794200247803e+01, 5.4644091485579487e+01, 1.7369227455294599e+01);
  mol.emplace_back(AtomicNumber(1), 8.2496252611433462e+01, 5.3255690904201302e+01, 2.0029980545066490e+01);
  mol.emplace_back(AtomicNumber(1), 7.4364327043788990e+01, 5.8968559336066981e+01, 1.9890858917756312e+01);
  mol.emplace_back(AtomicNumber(1), 7.6935280354563602e+01, 6.2670060117020732e+01, 1.6176073363099889e+01);
  mol.emplace_back(AtomicNumber(1), 7.3094695740819446e+01, 6.3539088407582163e+01, 1.4235797203894140e+01);
  mol.emplace_back(AtomicNumber(1), 7.3198082649677644e+01, 5.9929844049411386e+01, 1.5117165405163741e+01);
  mol.emplace_back(AtomicNumber(1), 6.9593657092778813e+01, 6.0008985773830709e+01, 1.6847511904251480e+01);
  mol.emplace_back(AtomicNumber(1), 7.0272371078988058e+01, 6.2954539467404786e+01, 1.8622814984617531e+01);
  mol.emplace_back(AtomicNumber(1), 6.6684197165594640e+01, 6.3637108494631597e+01, 1.5579864813570390e+01);
  mol.emplace_back(AtomicNumber(1), 6.9536757443250025e+01, 6.5488850991252690e+01, 1.4586114607734959e+01);
  mol.emplace_back(AtomicNumber(1), 7.0632023729214538e+01, 6.1677840589236382e+01, 1.1661801434277240e+01);
  mol.emplace_back(AtomicNumber(1), 6.7413196760371164e+01, 6.0484100681985090e+01, 1.2291911668049400e+01);
  mol.emplace_back(AtomicNumber(1), 6.5476511080544512e+01, 6.4094195416850908e+01, 1.0424578932019049e+01);
  mol.emplace_back(AtomicNumber(1), 6.8261721523951948e+01, 6.6087988616065132e+01, 1.0884765004860331e+01);
  mol.emplace_back(AtomicNumber(1), 7.4594174415831063e+01, 6.3126372251584556e+01, 2.2314829341106499e+01);
  mol.emplace_back(AtomicNumber(1), 7.5544423129399718e+01, 6.8720622583120701e+01, 2.1424617222208379e+01);
  mol.emplace_back(AtomicNumber(1), 7.1869302923252391e+01, 6.6018201035291369e+01, 2.5278486609655200e+01);
  mol.emplace_back(AtomicNumber(1), 7.2517422245699720e+01, 6.9548190285483486e+01, 2.5054157237501009e+01);
  mol.emplace_back(AtomicNumber(1), 7.0405351096833982e+01, 6.6494771032457280e+01, 2.0559557356223849e+01);
  mol.emplace_back(AtomicNumber(1), 6.8342677385320712e+01, 6.7801875601788694e+01, 2.3134327913496239e+01);
  mol.emplace_back(AtomicNumber(1), 6.9160399615280781e+01, 6.8771134958806684e+01, 1.7713762297609080e+01);
  mol.emplace_back(AtomicNumber(1), 7.7587349204327936e+01, 6.3346204075884927e+01, 2.4273303561586321e+01);
  mol.emplace_back(AtomicNumber(1), 8.1238469890414947e+01, 6.3070001725332695e+01, 2.6718722374911660e+01);
  mol.emplace_back(AtomicNumber(1), 8.0989026059866958e+01, 6.6683252302600152e+01, 2.7261243809093671e+01);
  mol.emplace_back(AtomicNumber(1), 7.5416488679944408e+01, 6.4125905018946341e+01, 2.8942892069444881e+01);
  mol.emplace_back(AtomicNumber(1), 7.7188314664490704e+01, 6.4943551659866856e+01, 3.4412004745989449e+01);
  mol.emplace_back(AtomicNumber(1), 7.1994553961803319e+01, 6.3099179094602853e+01, 3.2005230829139165e+01);
  mol.emplace_back(AtomicNumber(1), 7.3269627675621180e+01, 6.8607843736097180e+01, 3.3712031342403961e+01);
  mol.emplace_back(AtomicNumber(1), 7.3902175655919137e+01, 6.7454487273230825e+01, 3.0445375408779000e+01);
  mol.emplace_back(AtomicNumber(1), 7.1912237497722472e+01, 6.2023320294545378e+01, 3.6409350630062995e+01);
  mol.emplace_back(AtomicNumber(1), 7.3433259049008683e+01, 6.5050453459064585e+01, 3.7356802549167931e+01);
  mol.emplace_back(AtomicNumber(1), 7.0146666506199779e+01, 6.5050018822087111e+01, 3.6100966245918087e+01);
  mol.emplace_back(AtomicNumber(1), 6.8895535620662542e+01, 6.5627632467884851e+01, 3.1292671850466931e+01);
  mol.emplace_back(AtomicNumber(1), 6.8861709525459446e+01, 6.8649077557177179e+01, 3.3093939765921839e+01);
  mol.emplace_back(AtomicNumber(1), 6.9687085145674985e+01, 6.8648548433900245e+01, 2.9673762492950523e+01);
  mol.emplace_back(AtomicNumber(1), 7.9445252207673178e+01, 5.6798133448660920e+01, 3.5716690466054935e+01);
  mol.emplace_back(AtomicNumber(1), 7.8329161141309882e+01, 5.9177752003049278e+01, 4.1117168294679026e+01);
  mol.emplace_back(AtomicNumber(1), 8.1072457462281292e+01, 5.7313896362838690e+01, 4.0225746751147945e+01);
  mol.emplace_back(AtomicNumber(1), 8.1726888469531886e+01, 6.2310426364054145e+01, 4.0600611695585876e+01);
  mol.emplace_back(AtomicNumber(1), 8.2866242062819765e+01, 6.0946951268470855e+01, 3.7679870104247371e+01);
  mol.emplace_back(AtomicNumber(1), 7.9872065719548729e+01, 6.4225058941589168e+01, 3.6153160477734275e+01);
  mol.emplace_back(AtomicNumber(1), 7.7679170990133457e+01, 6.3809678271947071e+01, 3.8711925055879831e+01);
  mol.emplace_back(AtomicNumber(1), 7.1621975585812081e+01, 5.3843187816921507e+01, 3.6665389604312608e+01);
  mol.emplace_back(AtomicNumber(1), 7.4734429878734630e+01, 4.9035441442007162e+01, 3.8015183083735529e+01);
  mol.emplace_back(AtomicNumber(1), 7.2087282816083544e+01, 4.9310566648745670e+01, 3.5853355449579418e+01);
  mol.emplace_back(AtomicNumber(1), 7.6780114953606812e+01, 4.9048764010229611e+01, 3.3793724196908428e+01);
  mol.emplace_back(AtomicNumber(1), 7.4841860801209293e+01, 5.1567315219329252e+01, 3.2572696646375974e+01);
  mol.emplace_back(AtomicNumber(1), 7.8665513470091994e+01, 5.3878317823057017e+01, 3.4230364283926768e+01);
  mol.emplace_back(AtomicNumber(1), 7.9055911962159513e+01, 5.1841514460333151e+01, 3.6921088427884200e+01);
  mol.emplace_back(AtomicNumber(1), 7.6356381695093333e+01, 5.3186961569981371e+01, 4.1419297685800352e+01);
  mol.emplace_back(AtomicNumber(1), 7.2364505618669853e+01, 5.1500040974120850e+01, 4.5203946410270099e+01);
  mol.emplace_back(AtomicNumber(1), 7.6069861440641162e+01, 5.0288499848053164e+01, 4.6905569074325037e+01);
  mol.emplace_back(AtomicNumber(1), 7.7907733348502987e+01, 5.2590318383474397e+01, 4.4909772765562465e+01);
  mol.emplace_back(AtomicNumber(1), 7.6897088992325905e+01, 5.7028396560420568e+01, 4.6220618992352101e+01);
  mol.emplace_back(AtomicNumber(1), 7.3619661620563647e+01, 5.7258999822858243e+01, 4.2339272989025218e+01);
  mol.emplace_back(AtomicNumber(1), 7.1187036252183844e+01, 5.9777380956618870e+01, 4.6959690826650004e+01);
  mol.emplace_back(AtomicNumber(1), 7.3014571358885846e+01, 6.1942288844137046e+01, 4.1668080112252206e+01);
  mol.emplace_back(AtomicNumber(1), 7.1763799521286529e+01, 6.3908888886369574e+01, 4.4282591607073265e+01);
  mol.emplace_back(AtomicNumber(1), 7.5193236451603951e+01, 6.2172230702478565e+01, 4.7189859452110205e+01);
  mol.emplace_back(AtomicNumber(1), 7.6453532508187834e+01, 6.0282504713478573e+01, 4.4454122035094791e+01);
  mol.emplace_back(AtomicNumber(1), 7.3837849383253584e+01, 6.6035321952751701e+01, 4.2710925399281848e+01);
  mol.emplace_back(AtomicNumber(1), 7.6899583430631395e+01, 6.7790745115713477e+01, 4.2397892289204002e+01);
  mol.emplace_back(AtomicNumber(1), 6.9372294590427359e+01, 5.6842730982001321e+01, 4.1440859459334838e+01);
  mol.emplace_back(AtomicNumber(1), 6.5154955306256284e+01, 6.0440259039040285e+01, 4.0470296191384435e+01);
  mol.emplace_back(AtomicNumber(1), 6.6306044097935839e+01, 5.5036341909116217e+01, 3.8158235341102831e+01);
  mol.emplace_back(AtomicNumber(1), 6.4059027616195621e+01, 5.7500204448094202e+01, 3.6870727230277353e+01);
  mol.emplace_back(AtomicNumber(1), 6.7670369570214191e+01, 6.0420246840816780e+01, 3.6056368712577694e+01);
  mol.emplace_back(AtomicNumber(1), 6.9926891373679084e+01, 5.8033919759167468e+01, 3.7346465748008100e+01);
  mol.emplace_back(AtomicNumber(1), 7.2039075906104159e+01, 5.7730864402311539e+01, 3.3397807704953038e+01);
  mol.emplace_back(AtomicNumber(1), 7.0444260554947505e+01, 5.6023081231532458e+01, 3.0734371204276766e+01);
  mol.emplace_back(AtomicNumber(1), 6.1414771834307693e+01, 6.0387403403127955e+01, 4.1507302225108077e+01);
  mol.emplace_back(AtomicNumber(1), 5.9038724862038549e+01, 5.5211538405556411e+01, 4.3513227465171688e+01);
  mol.emplace_back(AtomicNumber(1), 5.9554563365255881e+01, 5.7520216646317706e+01, 4.7053326749404953e+01);
  mol.emplace_back(AtomicNumber(1), 5.8579767111090121e+01, 6.0569686372026908e+01, 4.5311811969722214e+01);
  mol.emplace_back(AtomicNumber(1), 5.4290315883178799e+01, 5.9330857602678066e+01, 4.5184179876425162e+01);
  mol.emplace_back(AtomicNumber(1), 5.5030937292787684e+01, 5.6113315647507207e+01, 4.6448387665806266e+01);
  mol.emplace_back(AtomicNumber(1), 5.6450593942023929e+01, 5.7912126919176416e+01, 5.0487847042632779e+01);
  mol.emplace_back(AtomicNumber(1), 5.6048384662525173e+01, 6.1199286379782023e+01, 4.9223884917630237e+01);
  mol.emplace_back(AtomicNumber(1), 5.1800016077614714e+01, 5.7381529655984998e+01, 5.0394967010273433e+01);
  mol.emplace_back(AtomicNumber(1), 5.4202934953447439e+01, 6.3868089702267056e+01, 4.9611996841251063e+01);
  mol.emplace_back(AtomicNumber(1), 5.0986659114689218e+01, 6.5223098825419626e+01, 5.0222869664455196e+01);
  mol.emplace_back(AtomicNumber(1), 4.7918122053751020e+01, 5.9178602379744333e+01, 5.0554629959084039e+01);
  mol.emplace_back(AtomicNumber(1), 5.6612089925043875e+01, 5.3760701277501660e+01, 4.0462151472371851e+01);
  mol.emplace_back(AtomicNumber(1), 5.4482822269678230e+01, 5.7902829467310539e+01, 3.6874090942537769e+01);
  mol.emplace_back(AtomicNumber(1), 5.5252337589658922e+01, 5.2163901714056550e+01, 3.6259570948174861e+01);
  mol.emplace_back(AtomicNumber(1), 5.3537297871082082e+01, 5.4498053461149574e+01, 3.3986797501204563e+01);
  mol.emplace_back(AtomicNumber(1), 5.9151182455643941e+01, 5.4963152821562247e+01, 3.5294072145874978e+01);
  mol.emplace_back(AtomicNumber(1), 5.6253287651512437e+01, 5.1406764099303807e+01, 3.1704518732509591e+01);
  mol.emplace_back(AtomicNumber(1), 5.9230777714300622e+01, 5.2946191581722879e+01, 3.0637069213103160e+01);
  mol.emplace_back(AtomicNumber(1), 5.9230248591023695e+01, 5.0921312389989602e+01, 3.3514328209434780e+01);
  mol.emplace_back(AtomicNumber(1), 5.4968746410489693e+01, 5.7367337813807609e+01, 3.2126494545853291e+01);
  mol.emplace_back(AtomicNumber(1), 5.7986336458764448e+01, 5.8700350526448211e+01, 3.3348013425142895e+01);
  mol.emplace_back(AtomicNumber(1), 5.7985901821786975e+01, 5.6726493936417931e+01, 3.0435511039116420e+01);
  mol.emplace_back(AtomicNumber(1), 5.0574150828550408e+01, 5.8695569519696036e+01, 3.6638857851427048e+01);
  mol.emplace_back(AtomicNumber(1), 4.7084639509042681e+01, 5.5188238084112044e+01, 3.9675439645891259e+01);
  mol.emplace_back(AtomicNumber(1), 4.6637322470186490e+01, 6.1114853422593505e+01, 3.8940959845746633e+01);
  mol.emplace_back(AtomicNumber(1), 5.0391414325414111e+01, 6.0278063857404419e+01, 4.1509097464797627e+01);
  mol.emplace_back(AtomicNumber(1), 4.7750578947566275e+01, 6.1503551161270927e+01, 4.3548546443906098e+01);
  mol.emplace_back(AtomicNumber(1), 4.3276010161772525e+01, 5.7086448942802654e+01, 4.1436853240238165e+01);
  mol.emplace_back(AtomicNumber(1), 4.3671567605789996e+01, 6.0124391339978942e+01, 4.3165782444834150e+01);
  mol.emplace_back(AtomicNumber(1), 4.2453922664777849e+01, 6.0123937805741583e+01, 3.9864846881768727e+01);
  mol.emplace_back(AtomicNumber(1), 4.6751065077464396e+01, 5.6807109647108668e+01, 4.4673349147078682e+01);
  mol.emplace_back(AtomicNumber(1), 4.9786380755515978e+01, 5.5952122020645511e+01, 4.3114306308893788e+01);
  mol.emplace_back(AtomicNumber(1), 4.9785870529498951e+01, 5.8077156689795792e+01, 4.5918414012191221e+01);
  mol.emplace_back(AtomicNumber(1), 4.3696228529946445e+01, 5.3851294741414321e+01, 3.7935871283977200e+01);
  mol.emplace_back(AtomicNumber(1), 4.1605738051875086e+01, 5.6792407578914244e+01, 3.3317078610702957e+01);
  mol.emplace_back(AtomicNumber(1), 4.4330817414312541e+01, 5.3680123361330701e+01, 3.1392846225143824e+01);
  mol.emplace_back(AtomicNumber(1), 4.2701136618658829e+01, 5.1117484844907686e+01, 3.3323144631127647e+01);
  mol.emplace_back(AtomicNumber(1), 4.1549235244803995e+01, 5.6112257400953368e+01, 2.8329562602455042e+01);
  mol.emplace_back(AtomicNumber(1), 3.9217729114095569e+01, 4.9121537358066000e+01, 3.1860440023861980e+01);
  mol.emplace_back(AtomicNumber(1), 3.8217723915236547e+01, 5.5580375124089429e+01, 2.4983632663591528e+01);
  mol.emplace_back(AtomicNumber(1), 3.5910368482667550e+01, 4.8506110295228375e+01, 2.8498636386690873e+01);
  mol.emplace_back(AtomicNumber(1), 3.5710643342890137e+01, 5.1537022911725579e+01, 2.4923766144260011e+01);
  mol.emplace_back(AtomicNumber(1), 3.7745424698805778e+01, 5.7088641024949894e+01, 3.3835921778242799e+01);
  mol.emplace_back(AtomicNumber(1), 3.3605677563743036e+01, 5.6792596551513149e+01, 3.5704407247126440e+01);
  mol.emplace_back(AtomicNumber(1), 3.1607292330375543e+01, 5.2972742231868331e+01, 3.4734580972311747e+01);
  mol.emplace_back(AtomicNumber(1), 3.4416785752741617e+01, 5.2916863034373606e+01, 3.2618371323530098e+01);
  mol.emplace_back(AtomicNumber(1), 3.4416483396583381e+01, 5.1037454949273538e+01, 3.5592705543916651e+01);
  mol.emplace_back(AtomicNumber(1), 3.7155943573797117e+01, 5.7753654497738879e+01, 3.9004549125276476e+01);
  mol.emplace_back(AtomicNumber(1), 3.7868351374390230e+01, 5.8238614878295948e+01, 4.3165858033873711e+01);
  mol.emplace_back(AtomicNumber(1), 3.4689076370496629e+01, 5.6543133823705261e+01, 4.3682679194605321e+01);
  mol.emplace_back(AtomicNumber(1), 3.9393322452993452e+01, 5.3125564372598760e+01, 4.0651464221949873e+01);
  mol.emplace_back(AtomicNumber(1), 4.0195983566821205e+01, 5.0400749571799764e+01, 4.5790290590137019e+01);
  mol.emplace_back(AtomicNumber(1), 4.0669246543506361e+01, 4.8505845733589908e+01, 4.0243302305585757e+01);
  mol.emplace_back(AtomicNumber(1), 4.1283709846089600e+01, 4.6670279291434653e+01, 4.3273402339907697e+01);
  mol.emplace_back(AtomicNumber(1), 3.6819232197077099e+01, 4.6831340637477119e+01, 4.4447300124274499e+01);
  mol.emplace_back(AtomicNumber(1), 3.6024885877600951e+01, 4.8984456634823943e+01, 4.1840253047109996e+01);
  mol.emplace_back(AtomicNumber(1), 3.9067552589749745e+01, 4.4536608574514638e+01, 4.0038380419338601e+01);
  mol.emplace_back(AtomicNumber(1), 3.6137759210923925e+01, 4.3539513553678681e+01, 4.1764909671928564e+01);
  mol.emplace_back(AtomicNumber(1), 3.6426036910545868e+01, 4.7367266927957523e+01, 3.7266341366074499e+01);
  mol.emplace_back(AtomicNumber(1), 3.5330581651982456e+01, 4.4071187960683829e+01, 3.7074779842569569e+01);
  mol.emplace_back(AtomicNumber(1), 3.1788158004782733e+01, 4.4883543368835149e+01, 3.9801956800854811e+01);
  mol.emplace_back(AtomicNumber(1), 3.3208211496476672e+01, 4.7749747468131119e+01, 4.1110629842757085e+01);
  mol.emplace_back(AtomicNumber(1), 4.4415722802998310e+01, 4.9741953297514591e+01, 4.6830924897759544e+01);
  mol.emplace_back(AtomicNumber(1), 4.8425551276317300e+01, 5.2538067459878555e+01, 4.3547752758990718e+01);
  mol.emplace_back(AtomicNumber(1), 4.8363379291279202e+01, 5.2665151532638795e+01, 4.8300149059687257e+01);
  mol.emplace_back(AtomicNumber(1), 4.9124957762106085e+01, 4.9232634348959316e+01, 4.8283784032622520e+01);
  mol.emplace_back(AtomicNumber(1), 5.2994341799922594e+01, 5.0057839893835833e+01, 4.5617626326522092e+01);
  mol.emplace_back(AtomicNumber(1), 5.2206288267989819e+01, 5.3480152557174719e+01, 4.5468375767910871e+01);
  mol.emplace_back(AtomicNumber(1), 5.0919611636599498e+01, 5.5117883585541570e+01, 4.9818430508289424e+01);
  mol.emplace_back(AtomicNumber(1), 5.3316558978306986e+01, 5.4469065064478308e+01, 5.2345787840497799e+01);
  mol.emplace_back(AtomicNumber(1), 5.0660775867886173e+01, 5.0396176434906394e+01, 4.0677920385795865e+01);
  mol.emplace_back(AtomicNumber(1), 4.9111691885663312e+01, 4.4877760807308810e+01, 3.9631371235827778e+01);
  mol.emplace_back(AtomicNumber(1), 5.2918072459006559e+01, 4.8388928386650477e+01, 3.6940155763113211e+01);
  mol.emplace_back(AtomicNumber(1), 5.1165294912429395e+01, 4.5430562350870979e+01, 3.5797948683581943e+01);
  mol.emplace_back(AtomicNumber(1), 4.8986006210134924e+01, 5.0694205120631580e+01, 3.6920691585426511e+01);
  mol.emplace_back(AtomicNumber(1), 4.9701210805191749e+01, 4.7396104146549654e+01, 3.2004909575721030e+01);
  mol.emplace_back(AtomicNumber(1), 4.7885184129762749e+01, 5.0405020352534912e+01, 3.2157561641112451e+01);
  mol.emplace_back(AtomicNumber(1), 5.1346141689576690e+01, 5.0404491229257992e+01, 3.2790468669348336e+01);
  mol.emplace_back(AtomicNumber(1), 4.6675305962565389e+01, 4.5460476713276854e+01, 3.6033559719890462e+01);
  mol.emplace_back(AtomicNumber(1), 4.5516091349133120e+01, 4.8206267472553741e+01, 3.7902158572333441e+01);
  mol.emplace_back(AtomicNumber(1), 4.5206006211598108e+01, 4.8205984013655389e+01, 3.4397510547653823e+01);
  mol.emplace_back(AtomicNumber(1), 5.1226200781054857e+01, 4.1463214917784704e+01, 4.0024112988121644e+01);
  mol.emplace_back(AtomicNumber(1), 5.6010703526304503e+01, 4.2147352417582375e+01, 4.3443213117279235e+01);
  mol.emplace_back(AtomicNumber(1), 5.2443203215230753e+01, 3.9729958240673909e+01, 4.5180343732667488e+01);
  mol.emplace_back(AtomicNumber(1), 5.2418202140396282e+01, 3.7750262397337615e+01, 4.2111674390910061e+01);
  mol.emplace_back(AtomicNumber(1), 5.5311561602154178e+01, 3.5039034726399542e+01, 4.3509485807713467e+01);
  mol.emplace_back(AtomicNumber(1), 5.7250269288789056e+01, 3.7575859585812807e+01, 4.5040390625922150e+01);
  mol.emplace_back(AtomicNumber(1), 5.1967219033121431e+01, 3.8961482269987165e+01, 4.6665101442224788e+01);
  mol.emplace_back(AtomicNumber(1), 5.9655871575526170e+01, 4.0477496047402532e+01, 4.2393961659146882e+01);
  mol.emplace_back(AtomicNumber(1), 6.0730804409849043e+01, 4.0480557403504704e+01, 3.7162935559956424e+01);
  mol.emplace_back(AtomicNumber(1), 6.3333089377521262e+01, 3.8333223967684233e+01, 4.2122294650968243e+01);
  mol.emplace_back(AtomicNumber(1), 6.4970310179871078e+01, 3.7993129981443900e+01, 3.8984895974990877e+01);
  mol.emplace_back(AtomicNumber(1), 6.7609841646786506e+01, 4.0631168564828009e+01, 4.2298341524103485e+01);
  mol.emplace_back(AtomicNumber(1), 5.9370541848447061e+01, 3.4851063682273711e+01, 4.0522377039641277e+01);
  mol.emplace_back(AtomicNumber(1), 6.1327485390875793e+01, 3.1743144731724751e+01, 3.6142086683438727e+01);
  mol.emplace_back(AtomicNumber(1), 5.9592414576815550e+01, 3.0703719848735194e+01, 3.9141629848518541e+01);
  mol.emplace_back(AtomicNumber(1), 5.5648424156953325e+01, 3.4900555605925618e+01, 3.7236502592708185e+01);
  mol.emplace_back(AtomicNumber(1), 5.2602110273645764e+01, 3.0864025304382064e+01, 3.4469320135235819e+01);
  mol.emplace_back(AtomicNumber(1), 5.1319345375052670e+01, 3.6484826286063658e+01, 3.5643898220958661e+01);
  mol.emplace_back(AtomicNumber(1), 4.9143401487758730e+01, 3.3853836380618525e+01, 3.4541866715953532e+01);
  mol.emplace_back(AtomicNumber(1), 5.2115449139698583e+01, 3.3622401638745693e+01, 3.9629349229019546e+01);
  mol.emplace_back(AtomicNumber(1), 4.8961987792814725e+01, 3.5032723041596277e+01, 3.9398065665225836e+01);
  mol.emplace_back(AtomicNumber(1), 4.7256018758985093e+01, 3.0868031523478738e+01, 3.9344737597816263e+01);
  mol.emplace_back(AtomicNumber(1), 4.8778136351344919e+01, 3.0312206418334171e+01, 3.6295419050186190e+01);
  mol.emplace_back(AtomicNumber(1), 5.1475455639003961e+01, 3.0037950485550599e+01, 4.1427725863711288e+01);
  mol.emplace_back(AtomicNumber(1), 5.0654501977602685e+01, 2.6597798808875549e+01, 3.5454717752199869e+01);
  mol.emplace_back(AtomicNumber(1), 5.3008835998258228e+01, 2.4235811397964561e+01, 3.6651348940214341e+01);
  mol.emplace_back(AtomicNumber(1), 5.4239860199272499e+01, 2.7432490778216849e+01, 4.2451447123732258e+01);
  mol.emplace_back(AtomicNumber(1), 5.0450940694067604e+01, 3.2169693681961832e+01, 3.0644023404742679e+01);
  mol.emplace_back(AtomicNumber(1), 5.4073658798539952e+01, 3.5155404052802162e+01, 2.7107892545086319e+01);
  mol.emplace_back(AtomicNumber(1), 5.1856859035363826e+01, 3.2987756062599928e+01, 2.3575371062068950e+01);
  mol.emplace_back(AtomicNumber(1), 4.8914215519812807e+01, 2.9313070493430089e+01, 2.6313035999593140e+01);
  mol.emplace_back(AtomicNumber(1), 5.2206741802227178e+01, 2.8158334530591752e+01, 2.5659115218359581e+01);
  mol.emplace_back(AtomicNumber(1), 5.4600438815233588e+01, 3.0005409404020021e+01, 2.3860984248046410e+01);
  mol.emplace_back(AtomicNumber(1), 5.4599853000176999e+01, 3.0079222101150361e+01, 2.7378557998490791e+01);
  mol.emplace_back(AtomicNumber(1), 5.1832689439964518e+01, 3.7668172700375465e+01, 2.4022782587224590e+01);
  mol.emplace_back(AtomicNumber(1), 4.8138558590367872e+01, 4.1303873222392234e+01, 2.6312166725638200e+01);
  mol.emplace_back(AtomicNumber(1), 4.9419849502689537e+01, 3.9952227811500087e+01, 2.0737512852607981e+01);
  mol.emplace_back(AtomicNumber(1), 4.7850904500322287e+01, 4.2948936387596397e+01, 2.2142146180231681e+01);
  mol.emplace_back(AtomicNumber(1), 5.3316899128985007e+01, 4.1247710565999149e+01, 2.2940895561262199e+01);
  mol.emplace_back(AtomicNumber(1), 5.0522240055632587e+01, 4.4632436579416826e+01, 1.9031978455755809e+01);
  mol.emplace_back(AtomicNumber(1), 5.3512995994863537e+01, 4.5896757752357281e+01, 2.0385362414557832e+01);
  mol.emplace_back(AtomicNumber(1), 5.3512447974326733e+01, 4.2810646239721379e+01, 1.8695758407792933e+01);
  mol.emplace_back(AtomicNumber(1), 4.9689229942421491e+01, 4.5302268853477770e+01, 2.4935066705674227e+01);
  mol.emplace_back(AtomicNumber(1), 5.2734334401096085e+01, 4.4154317006939941e+01, 2.6270876212778550e+01);
  mol.emplace_back(AtomicNumber(1), 5.2733861969598841e+01, 4.6581254300092866e+01, 2.3723563374126329e+01);
  mol.emplace_back(AtomicNumber(1), 4.7548302677703724e+01, 3.5453300457708117e+01, 2.2544109795351872e+01);
  mol.emplace_back(AtomicNumber(1), 4.2099032124043653e+01, 3.6551004490198437e+01, 2.1678483011570638e+01);
  mol.emplace_back(AtomicNumber(1), 4.1631249352726591e+01, 3.1850178811741706e+01, 2.0817637234541582e+01);
  mol.emplace_back(AtomicNumber(1), 4.4401379782741806e+01, 3.3307762261577196e+01, 1.9170759932387970e+01);
  mol.emplace_back(AtomicNumber(1), 4.4063421186869043e+01, 2.9514855434535509e+01, 2.3160859666381800e+01);
  mol.emplace_back(AtomicNumber(1), 4.5491222555117879e+01, 3.4030469068810348e+01, 2.6929861459662629e+01);
  mol.emplace_back(AtomicNumber(1), 4.1247918435857940e+01, 3.1507647078975570e+01, 2.9389018580927999e+01);
  mol.emplace_back(AtomicNumber(1), 4.5813572014321501e+01, 3.4315647617810342e+01, 3.1559916899831311e+01);
  mol.emplace_back(AtomicNumber(1), 4.3245774545948521e+01, 3.2622094083728427e+01, 3.3490158614035472e+01);
  mol.emplace_back(AtomicNumber(1), 4.8625503183213389e+01, 3.1578341728224064e+01, 3.3413738095040308e+01);
  mol.emplace_back(AtomicNumber(1), 4.3076681864452802e+01, 3.7962459728642429e+01, 2.9106258881193927e+01);
  mol.emplace_back(AtomicNumber(1), 3.8568060216337258e+01, 3.9524828484567962e+01, 3.2396668670500617e+01);
  mol.emplace_back(AtomicNumber(1), 4.3023089235404754e+01, 4.2332186516566473e+01, 2.9629618493847481e+01);
  mol.emplace_back(AtomicNumber(1), 4.0597285777845244e+01, 4.3823010246548456e+01, 3.1800346737411779e+01);
  mol.emplace_back(AtomicNumber(1), 4.6688760811607068e+01, 4.0489835958110703e+01, 3.1442659402213856e+01);
  mol.emplace_back(AtomicNumber(1), 4.0451266650675208e+01, 4.2567362915897519e+01, 3.6252767934614461e+01);
  mol.emplace_back(AtomicNumber(1), 4.9303291203688019e+01, 3.9560298641381486e+01, 3.5269600194317427e+01);
  mol.emplace_back(AtomicNumber(1), 4.3052512269053487e+01, 4.1671179262874162e+01, 4.0060905953127481e+01);
  mol.emplace_back(AtomicNumber(1), 4.9645652861115153e+01, 3.9716994720389373e+01, 3.9946048407516059e+01);
  mol.emplace_back(AtomicNumber(1), 3.8829711676774195e+01, 3.7727188843011930e+01, 2.6293307260267980e+01);
  mol.emplace_back(AtomicNumber(1), 3.6217373366840491e+01, 3.9005626269090207e+01, 2.2891479226649849e+01);
  mol.emplace_back(AtomicNumber(1), 3.2972770535507166e+01, 3.8426368561682040e+01, 2.6606548240204621e+01);
  mol.emplace_back(AtomicNumber(1), 3.2848558846250192e+01, 4.2010441769679332e+01, 2.6567846651949903e+01);
  mol.emplace_back(AtomicNumber(1), 3.3797409162586980e+01, 3.7964273865591871e+01, 2.1101587461648720e+01);
  mol.emplace_back(AtomicNumber(1), 3.0699543554299499e+01, 3.8809321533352886e+01, 1.9603922923586552e+01);
  mol.emplace_back(AtomicNumber(1), 4.0516292121956702e+01, 4.1927369415202890e+01, 2.3568568048508549e+01);
  mol.emplace_back(AtomicNumber(1), 3.8835513135560433e+01, 4.7527099952107136e+01, 2.2721951908176660e+01);
  mol.emplace_back(AtomicNumber(1), 4.3891021484892569e+01, 4.4518297129681230e+01, 2.2578975239848919e+01);
  mol.emplace_back(AtomicNumber(1), 4.1737697617686962e+01, 4.8722710688087552e+01, 2.6262240165008819e+01);
  mol.emplace_back(AtomicNumber(1), 4.1896434600762959e+01, 4.5267441203500496e+01, 2.6940859664918609e+01);
  mol.emplace_back(AtomicNumber(1), 4.2720147262108171e+01, 4.9778991926898989e+01, 2.0669841764941889e+01);
  mol.emplace_back(AtomicNumber(1), 4.5408206892421113e+01, 4.9576129841979835e+01, 2.2929935150525999e+01);
  mol.emplace_back(AtomicNumber(1), 4.5407621077364524e+01, 4.7605996012147891e+01, 2.0014900531674272e+01);
  mol.emplace_back(AtomicNumber(1), 4.6463883418916069e+01, 4.5262716888528004e+01, 2.7122745791359861e+01);
  mol.emplace_back(AtomicNumber(1), 4.6726914379324981e+01, 4.8292363388612578e+01, 2.5354472491672890e+01);
  mol.emplace_back(AtomicNumber(1), 4.5595573224230456e+01, 4.8291853162595551e+01, 2.8685964923980439e+01);
  mol.emplace_back(AtomicNumber(1), 3.7288866899863379e+01, 4.8423548166768960e+01, 1.8996924038659859e+01);
  mol.emplace_back(AtomicNumber(1), 3.8772301801228380e+01, 4.4975932072437359e+01, 1.4417532152256269e+01);
  mol.emplace_back(AtomicNumber(1), 3.4210559955562047e+01, 4.8205984013655389e+01, 1.5877591143137341e+01);
  mol.emplace_back(AtomicNumber(1), 3.5068023123070802e+01, 4.7093861371868996e+01, 1.2477652835508211e+01);
  mol.emplace_back(AtomicNumber(1), 3.4843372497498478e+01, 4.2636545784354809e+01, 1.4182072294026870e+01);
  mol.emplace_back(AtomicNumber(1), 3.3249124064138520e+01, 4.3962320846457544e+01, 1.7016396715888408e+01);
  mol.emplace_back(AtomicNumber(1), 3.2983825432542808e+01, 4.6855566924656095e+01, 1.0813749102193709e+01);
  mol.emplace_back(AtomicNumber(1), 2.9570148816973649e+01, 4.6120312336855974e+01, 1.0213798895205990e+01);
  mol.emplace_back(AtomicNumber(1), 3.8625224427504506e+01, 4.7337125798432972e+01, 1.0354470097827150e+01);
  mol.emplace_back(AtomicNumber(1), 4.2987033263534641e+01, 5.0503739638200265e+01, 1.0478095972027530e+01);
  mol.emplace_back(AtomicNumber(1), 4.1721257001582664e+01, 5.2089068564892152e+01, 6.1073865265892096e+00);
  mol.emplace_back(AtomicNumber(1), 4.2015430646290291e+01, 4.8567809848769222e+01, 6.6240376119818096e+00);
  mol.emplace_back(AtomicNumber(1), 3.7510059326875833e+01, 4.7949075765450843e+01, 6.4936843132605899e+00);
  mol.emplace_back(AtomicNumber(1), 3.6822577012077630e+01, 5.1357291072911785e+01, 6.8581179702392401e+00);
  mol.emplace_back(AtomicNumber(1), 3.8120214054204155e+01, 5.2345296511740663e+01, 2.6811243359333101e+00);
  mol.emplace_back(AtomicNumber(1), 4.0177539841168560e+01, 4.9390974489577616e+01, 2.5060790175722398e+00);
  mol.emplace_back(AtomicNumber(1), 3.6795610622214596e+01, 4.8840402822682471e+01, -4.0940913551685004e-01);
  mol.emplace_back(AtomicNumber(1), 3.6144467738184872e+01, 4.6842036486574855e+01, 2.4171863070496800e+00);
  mol.emplace_back(AtomicNumber(1), 3.2267713769011259e+01, 4.8712166017068931e+01, 2.2447866050732097e+00);
  mol.emplace_back(AtomicNumber(1), 3.3867499099518987e+01, 5.1371104969891377e+01, 3.7663183823764501e+00);
  mol.emplace_back(AtomicNumber(1), 4.3670263694857589e+01, 5.4478475899903529e+01, 1.0216954737607621e+01);
  mol.emplace_back(AtomicNumber(1), 4.3400713179786628e+01, 5.8547547282977675e+01, 1.1050286104236839e+01);
  mol.emplace_back(AtomicNumber(1), 3.7987725365115359e+01, 5.7784740490257931e+01, 8.6806264058105089e+00);
  mol.emplace_back(AtomicNumber(1), 3.9090588349555645e+01, 6.0956891227172996e+01, 9.6861873990771894e+00);
  mol.emplace_back(AtomicNumber(1), 4.1660710180895101e+01, 6.1311876254206652e+01, 6.1460314230642599e+00);
  mol.emplace_back(AtomicNumber(1), 4.3169713074891270e+01, 5.8178257030207291e+01, 6.8576644360018797e+00);
  mol.emplace_back(AtomicNumber(1), 4.3455212877309386e+01, 5.6973140972502215e+01, 3.3073039396483503e+00);
  mol.emplace_back(AtomicNumber(1), 4.1673617009399969e+01, 5.4931537705766281e+01, 1.4892533676851309e+01);
  mol.emplace_back(AtomicNumber(1), 3.8483929615306984e+01, 5.7925884124376338e+01, 1.8496808055671011e+01);
  mol.emplace_back(AtomicNumber(1), 4.2731617898861401e+01, 5.3759850900806612e+01, 1.8902532225509308e+01);
  mol.emplace_back(AtomicNumber(1), 4.0995167584829190e+01, 5.5002591402952682e+01, 2.1687950538775532e+01);
  mol.emplace_back(AtomicNumber(1), 3.9450770123059044e+01, 5.1026872483735140e+01, 1.9177733021287381e+01);
  mol.emplace_back(AtomicNumber(1), 4.0188840402582777e+01, 5.9992904205664324e+01, 2.2102726496101141e+01);
  mol.emplace_back(AtomicNumber(1), 4.5650526455990580e+01, 6.1939983378430476e+01, 2.1266730615827427e+01);
  mol.emplace_back(AtomicNumber(1), 4.1067865343626018e+01, 6.4955419139077762e+01, 2.3772129332043630e+01);
  mol.emplace_back(AtomicNumber(1), 4.0472696143390472e+01, 6.5675177973768086e+01, 1.9612369998757380e+01);
  mol.emplace_back(AtomicNumber(1), 4.6251024683515105e+01, 6.5116707252238925e+01, 2.5054459593659249e+01);
  mol.emplace_back(AtomicNumber(1), 4.5878540793823319e+01, 6.7631592392919899e+01, 2.2623119238951961e+01);
  mol.emplace_back(AtomicNumber(1), 4.3828528246436335e+01, 6.7630987680603425e+01, 2.5482539221947420e+01);
  mol.emplace_back(AtomicNumber(1), 4.8311071675903683e+01, 6.0775968860986140e+01, 2.4331148074109610e+01);
  mol.emplace_back(AtomicNumber(1), 4.5813760986920400e+01, 5.8853834071274790e+01, 2.9278545199611060e+01);
  mol.emplace_back(AtomicNumber(1), 5.1530975788560781e+01, 5.8931426220383130e+01, 2.7954036253920961e+01);
  mol.emplace_back(AtomicNumber(1), 4.9544136780986072e+01, 5.6803349092390562e+01, 3.0177865797776160e+01);
  mol.emplace_back(AtomicNumber(1), 5.0058973729429226e+01, 5.6582572405095689e+01, 2.4447933140229811e+01);
  mol.emplace_back(AtomicNumber(1), 5.0949639382564712e+01, 5.1763619955066574e+01, 2.6689998539878861e+01);
  mol.emplace_back(AtomicNumber(1), 5.3362592703399031e+01, 5.4012847416213930e+01, 2.5467950537312340e+01);
  mol.emplace_back(AtomicNumber(1), 5.2528410960074758e+01, 5.4012715135394700e+01, 2.8885973522656201e+01);
  mol.emplace_back(AtomicNumber(1), 4.7397219084883169e+01, 5.2299243889388734e+01, 2.6367025471098870e+01);
  mol.emplace_back(AtomicNumber(1), 4.6027186640118060e+01, 5.5025457087419582e+01, 2.8117894394427150e+01);
  mol.emplace_back(AtomicNumber(1), 4.6048578338313540e+01, 5.5025173628521230e+01, 2.4599602548106951e+01);
  mol.emplace_back(AtomicNumber(1), 4.7532353390356562e+01, 6.0749890642337945e+01, 3.3094393300159197e+01);
  mol.emplace_back(AtomicNumber(1), 4.9771073975005081e+01, 6.6132416074066526e+01, 3.1991416932159570e+01);
  mol.emplace_back(AtomicNumber(1), 4.7132733035462728e+01, 6.3623880412708587e+01, 3.6715656315620009e+01);
  mol.emplace_back(AtomicNumber(1), 4.8573951458233466e+01, 6.6800528697477475e+01, 3.6748405267009382e+01);
  mol.emplace_back(AtomicNumber(1), 4.6485142836292319e+01, 7.0292799016929152e+01, 3.3926949879132927e+01);
  mol.emplace_back(AtomicNumber(1), 4.2632520667998243e+01, 6.3228889886487813e+01, 3.3937286680292758e+01);
  mol.emplace_back(AtomicNumber(1), 4.2118930938707820e+01, 7.1021760817185907e+01, 3.1770243402407012e+01);
  mol.emplace_back(AtomicNumber(1), 5.3734018346536978e+01, 6.6804723889173061e+01, 3.3261482872106583e+01);
  mol.emplace_back(AtomicNumber(1), 5.6290156205557828e+01, 6.2509792455893646e+01, 3.6343550371126021e+01);
  mol.emplace_back(AtomicNumber(1), 5.8190861502553922e+01, 6.3290759515367668e+01, 3.2056801451378973e+01);
  mol.emplace_back(AtomicNumber(1), 5.8256510583411774e+01, 6.6794840622250589e+01, 3.2695245376762621e+01);
  mol.emplace_back(AtomicNumber(1), 6.1370892396843125e+01, 6.5974434981386139e+01, 3.6328942789231050e+01);
  mol.emplace_back(AtomicNumber(1), 5.9609478802496220e+01, 6.1016077445148476e+01, 3.6000999741099989e+01);
  mol.emplace_back(AtomicNumber(1), 6.2625953912437474e+01, 6.0812232702715050e+01, 3.4202642003668139e+01);
  mol.emplace_back(AtomicNumber(1), 6.2625424789160547e+01, 6.2020920342539334e+01, 3.7506865689954424e+01);
  mol.emplace_back(AtomicNumber(1), 6.3792292792848272e+01, 6.3332938199442154e+01, 3.1767654477802079e+01);
  mol.emplace_back(AtomicNumber(1), 6.2081089218029099e+01, 6.6317439145649516e+01, 3.1033231369437122e+01);
  mol.emplace_back(AtomicNumber(1), 6.4679424658384320e+01, 6.6316891125112704e+01, 3.3405442197948595e+01);
  mol.emplace_back(AtomicNumber(1), 5.6545155830513487e+01, 6.3264945858357926e+01, 4.0305153037205727e+01);
  mol.emplace_back(AtomicNumber(1), 5.7541249296575280e+01, 6.8957424146802296e+01, 4.1900289641780518e+01);
  mol.emplace_back(AtomicNumber(1), 5.5825850530060528e+01, 6.4389559588931604e+01, 4.4898132053470235e+01);
  mol.emplace_back(AtomicNumber(1), 5.4541044727399317e+01, 6.8557047901512874e+01, 4.7792058433024827e+01);
  mol.emplace_back(AtomicNumber(1), 5.7506478338377683e+01, 6.6692568651725907e+01, 4.8115617316861410e+01);
  mol.emplace_back(AtomicNumber(1), 5.7505930317840871e+01, 6.9692773220901870e+01, 4.6277820998039132e+01);
  mol.emplace_back(AtomicNumber(1), 5.1975968464450503e+01, 6.5131031375235537e+01, 4.3289105962876285e+01);
  mol.emplace_back(AtomicNumber(1), 5.1998966429736633e+01, 6.8127078547235698e+01, 4.5132458076106232e+01);
  mol.emplace_back(AtomicNumber(1), 5.3010423368088986e+01, 6.8126530526698886e+01, 4.1762623103481872e+01);
  mol.emplace_back(AtomicNumber(1), 6.1337916678335070e+01, 6.9673724782932752e+01, 4.3658150551268101e+01);
  mol.emplace_back(AtomicNumber(1), 6.5088833793901173e+01, 6.5118559183708143e+01, 4.3560206053258227e+01);
  mol.emplace_back(AtomicNumber(1), 6.5204409435388399e+01, 7.0777173582429640e+01, 4.2588924689432005e+01);
  mol.emplace_back(AtomicNumber(1), 6.7705839727027708e+01, 6.9734781829637342e+01, 4.4840476513545845e+01);
  mol.emplace_back(AtomicNumber(1), 6.9534905511780806e+01, 6.6801265690613192e+01, 4.1799510554787155e+01);
  mol.emplace_back(AtomicNumber(1), 6.6603770427502795e+01, 6.4233392633200651e+01, 3.9938187147401820e+01);
  mol.emplace_back(AtomicNumber(1), 6.4098409505806387e+01, 6.6699296076246753e+01, 4.0069617589936769e+01);
  mol.emplace_back(AtomicNumber(1), 6.6430803807729632e+01, 6.6699107103647847e+01, 3.7435471842089996e+01);
  mol.emplace_back(AtomicNumber(1), 6.7468660218148315e+01, 7.1271129058694342e+01, 3.8649034974965915e+01);
  mol.emplace_back(AtomicNumber(1), 7.0256232819041998e+01, 7.1329672769833564e+01, 4.0793987356040255e+01);
  mol.emplace_back(AtomicNumber(1), 7.0255665901245308e+01, 6.9236593367157269e+01, 3.7965955721722082e+01);
  mol.emplace_back(AtomicNumber(1), 6.6836981511805277e+01, 6.3647445295791414e+01, 4.6779486556338959e+01);
  mol.emplace_back(AtomicNumber(1), 6.5350239589959543e+01, 6.5992368481021742e+01, 5.2028975278441948e+01);
  mol.emplace_back(AtomicNumber(1), 6.3750227492333124e+01, 6.2036718451807381e+01, 5.2650147108286141e+01);
  mol.emplace_back(AtomicNumber(1), 6.5372349384030827e+01, 6.0761947094147764e+01, 4.9753518420567275e+01);
  mol.emplace_back(AtomicNumber(1), 6.9537551128165418e+01, 6.1150077915028461e+01, 5.2579509150817316e+01);
  mol.emplace_back(AtomicNumber(1), 6.7245011147350169e+01, 6.1143841819264772e+01, 5.5267795548248941e+01);
  mol.emplace_back(AtomicNumber(1), 6.9232190305602899e+01, 5.6855883474884756e+01, 5.3189701672665421e+01);
  mol.emplace_back(AtomicNumber(1), 6.6002402926023322e+01, 5.7049467005197918e+01, 5.4553800377825070e+01);
  mol.emplace_back(AtomicNumber(1), 6.6470488053498627e+01, 5.8391418121766492e+01, 4.8892804924357435e+01);
  mol.emplace_back(AtomicNumber(1), 6.5163629148545795e+01, 5.4456951920888820e+01, 5.4254581164726808e+01);
  mol.emplace_back(AtomicNumber(1), 6.3400458109029124e+01, 5.1988119402559768e+01, 5.2424305955340749e+01);
  mol.emplace_back(AtomicNumber(1), 6.4836971114087262e+01, 5.5196515083943858e+01, 4.6438202042725557e+01);
  mol.emplace_back(AtomicNumber(1), 6.9096659157671823e+01, 6.8490756313818750e+01, 5.2433660098986302e+01);
  mol.emplace_back(AtomicNumber(1), 7.4155890267202295e+01, 6.5891343729649805e+01, 5.0882308445576641e+01);
  mol.emplace_back(AtomicNumber(1), 7.2174115725278099e+01, 7.1053413727501663e+01, 5.2319520649250698e+01);
  mol.emplace_back(AtomicNumber(1), 7.5648358058794713e+01, 7.0053427425902527e+01, 5.2951123766554169e+01);
  mol.emplace_back(AtomicNumber(1), 7.6697911873085317e+01, 6.9959148996311299e+01, 4.8863268507149371e+01);
  mol.emplace_back(AtomicNumber(1), 7.3274370887853564e+01, 6.7117435745832779e+01, 4.7324804784984686e+01);
  mol.emplace_back(AtomicNumber(1), 7.1034837721029788e+01, 6.9829910635923596e+01, 4.7282285950232193e+01);
  mol.emplace_back(AtomicNumber(1), 7.3838888732547545e+01, 6.9829646074285151e+01, 4.5157175692042344e+01);
  mol.emplace_back(AtomicNumber(1), 7.2239708114356290e+01, 7.3917142285752035e+01, 4.8826437747623757e+01);
  mol.emplace_back(AtomicNumber(1), 7.5268504237745816e+01, 7.4405107330631608e+01, 5.0547789151003862e+01);
  mol.emplace_back(AtomicNumber(1), 7.5268050703508464e+01, 7.4080754761879646e+01, 4.7044426139996759e+01);
  mol.emplace_back(AtomicNumber(1), 7.7259141594558429e+01, 6.5076739547571577e+01, 5.3172391782606176e+01);
  mol.emplace_back(AtomicNumber(1), 7.5644748682155722e+01, 6.3297562528928069e+01, 5.8498357715483998e+01);
  mol.emplace_back(AtomicNumber(1), 7.7898228026778327e+01, 6.0648733610146770e+01, 5.4727560682513619e+01);
  mol.emplace_back(AtomicNumber(1), 8.0774126420397877e+01, 6.2062021882800096e+01, 5.6382866162578168e+01);
  mol.emplace_back(AtomicNumber(1), 8.0091595187690842e+01, 5.9887476392738009e+01, 6.0070685327371557e+01);
  mol.emplace_back(AtomicNumber(1), 7.6625081833469253e+01, 5.9385243916641478e+01, 5.9353912259743865e+01);
  mol.emplace_back(AtomicNumber(1), 7.8020190942108385e+01, 5.5404139278355288e+01, 5.8668943280511030e+01);
  mol.emplace_back(AtomicNumber(1), 7.8640625778816869e+01, 5.6821414872845402e+01, 5.5465876626415920e+01);
  mol.emplace_back(AtomicNumber(1), 8.2536957309236513e+01, 5.5416762647961811e+01, 5.9518545187905538e+01);
  mol.emplace_back(AtomicNumber(1), 8.0696609860329090e+01, 5.7654784033994396e+01, 5.3270166205277043e+01);
  mol.emplace_back(AtomicNumber(1), 8.4025059833794344e+01, 5.7532480967986324e+01, 5.2061705332571428e+01);
  mol.emplace_back(AtomicNumber(1), 8.6613549801746885e+01, 5.5641696732432479e+01, 5.8100343627680822e+01);
  mol.emplace_back(AtomicNumber(1), 7.5814767329386058e+01, 6.6776415793857836e+01, 6.0660260938679663e+01);
  mol.emplace_back(AtomicNumber(1), 7.9018627668396547e+01, 6.8463487567797486e+01, 6.3834093531724953e+01);
  mol.emplace_back(AtomicNumber(1), 8.1071777160925251e+01, 6.9211856853961265e+01, 6.0970138206355891e+01);
  mol.emplace_back(AtomicNumber(1), 7.7726508626157909e+01, 7.1577472538771133e+01, 6.5464416834214916e+01);
  mol.emplace_back(AtomicNumber(1), 7.7725431482344163e+01, 7.6728733303965896e+01, 6.3542830065040377e+01);
  mol.emplace_back(AtomicNumber(1), 7.4517470437937547e+01, 7.5454755631221659e+01, 6.2946016803194404e+01);
  mol.emplace_back(AtomicNumber(1), 7.8916412389651526e+01, 7.3924266552730558e+01, 6.7862649189594933e+01);
#endif
  return mol;
}



BasisSet<double> make_631Gd( const Molecule& mol, SphericalType sph ) {

#if 1
  std::string basis_path = GAUXC_REF_DATA_PATH  "/../basis/old/6-31g*.g94";
  return parse_basis( mol, basis_path, sph );
#else
  BasisSet<double> basis;

  for( const auto& atom : mol ) {

    if( atom.Z == 1ll ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          {0.1873113696e+02, 0.2825394365e+01, 0.6401216923e+00},
          {0.3349460434e-01, 0.2347269535e+00, 0.8137573261e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {0.1612777588},
          {1.0000000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z.get() == 6 ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.3047524880e+04, 0.4573695180e+03, 0.1039486850e+03, 
            0.2921015530e+02, 0.9286662960e+01, 0.3163926960e+01 }, 
          { 0.1834737132e-02, 0.1403732281e-01, 0.6884262226e-01,
            0.2321844432e+00, 0.4679413484e+00, 0.3623119853e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 0.7868272350e+01, 0.1881288540e+01, 0.5442492580e+00},
          {-0.1193324198e+00, -0.1608541517e+00, 0.1143456438e+01},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 0.7868272350e+01, 0.1881288540e+01, 0.5442492580e+00},
          { 0.6899906659e-01, 0.3164239610e+00, 0.7443082909e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          { 0.1687144782 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          { 0.1687144782 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(2), sph,
          { 0.8 },
          { 1.0 },
          { atom.x, atom.y, atom.z }
        )
      );


    } else if( atom.Z.get() == 7 ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.4173511460e+04, 0.6274579110e+03, 0.1429020930e+03, 
            0.4023432930e+02, 0.1282021290e+02, 0.4390437010e+01 }, 
          { 0.1834772160e-02, 0.1399462700e-01, 0.6858655181e-01,
            0.2322408730e+00, 0.4690699481e+00, 0.3604551991e+00 },
          { atom.x, atom.y, atom.z }
        )
      );
      
      
      

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 0.1162636186e+02, 0.2716279807e+01, 0.7722183966e+00},
          {-0.1149611817e+00, -0.1691174786e+00, 0.1145851947e+01},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 0.1162636186e+02, 0.2716279807e+01, 0.7722183966e+00},
          { 0.6757974388e-01, 0.3239072959e+00, 0.7408951398e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          { 0.2120314975 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          { 0.2120314975 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(2), sph,
          { 0.8 },
          { 1.0 },
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z.get() == 8 ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.5484671660e+04, 0.8252349460e+03, 0.1880469580e+03, 
            0.5296450000e+02, 0.1689757040e+02, 0.5799635340e+01 }, 
          { 0.1831074430e-02, 0.1395017220e-01, 0.6844507810e-01,
            0.2327143360e+00, 0.4701928980e+00, 0.3585208530e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 0.1553961625e+02,  0.3599933586e+01, 0.1013761750e+01},
          {-0.1107775495e+00, -0.1480262627e+00, 0.1130767015e+01},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 0.1553961625e+02,  0.3599933586e+01, 0.1013761750e+01},
          { 0.7087426823e-01,  0.3397528391e+00, 0.7271585773e+00},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          { 0.2700058226 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          { 0.2700058226 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(2), sph,
          { 0.8 },
          { 1.0 },
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z.get() == 16 ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.2191710000e+05, 0.3301490000e+04, 0.7541460000e+03, 
            0.2127110000e+03, 0.6798960000e+02, 0.2305150000e+02 }, 
          { 0.1869240849e-02, 0.1423030646e-01, 0.6969623166e-01,
            0.2384871083e+00, 0.4833072195e+00, 0.3380741536e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(0), SphericalType(false),
          { 0.4237350000e+03, 0.1007100000e+03, 0.3215990000e+02, 
            0.1180790000e+02, 0.4631100000e+01, 0.1870250000e+01 }, 
          { -0.2376770499e-02, -0.3169300665e-01, -0.1133170238e+00,
             0.5609001177e-01,  0.5922551243e+00,  0.4550060955e+00 },
          { atom.x, atom.y, atom.z }
        )
      );


      basis.emplace_back(
        Shell<double>( 
          PrimSize(6), AngularMomentum(1), SphericalType(false),
          { 0.4237350000e+03, 0.1007100000e+03, 0.3215990000e+02, 
            0.1180790000e+02, 0.4631100000e+01, 0.1870250000e+01 }, 
          { 0.4061009982e-02, 0.3068129986e-01, 0.1304519994e+00,
            0.3272049985e+00, 0.4528509980e+00, 0.2560419989e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      
      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 0.2615840000e+01, 0.9221670000e+00, 0.3412870000e+00 },
          {-0.2503731142e+00, 0.6695676310e-01, 0.1054506269e+01 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 0.2615840000e+01, 0.9221670000e+00, 0.3412870000e+00 },
          {-0.1451048955e-01, 0.3102627765e+00, 0.7544824565e+00 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          { 0.1171670000 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          { 0.1171670000 },
          { 1.0000000000 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>( 
          PrimSize(1), AngularMomentum(2), sph,
          { 0.65 },
          { 1.0 },
          { atom.x, atom.y, atom.z }
        )
      );
      
    } else { throw std::runtime_error("NYI"); }

  }

  return basis;
#endif

}



/*
H     0
S   4   1.00
      1.301000D+01           1.968500D-02
      1.962000D+00           1.379770D-01
      4.446000D-01           4.781480D-01
      1.220000D-01           5.012400D-01
S   1   1.00
      1.220000D-01           1.000000D+00
P   1   1.00
      7.270000D-01           1.0000000
****
C     0
S   9   1.00
      6.665000D+03           6.920000D-04
      1.000000D+03           5.329000D-03
      2.280000D+02           2.707700D-02
      6.471000D+01           1.017180D-01
      2.106000D+01           2.747400D-01
      7.495000D+00           4.485640D-01
      2.797000D+00           2.850740D-01
      5.215000D-01           1.520400D-02
      1.596000D-01          -3.191000D-03
S   9   1.00
      6.665000D+03          -1.460000D-04
      1.000000D+03          -1.154000D-03
      2.280000D+02          -5.725000D-03
      6.471000D+01          -2.331200D-02
      2.106000D+01          -6.395500D-02
      7.495000D+00          -1.499810D-01
      2.797000D+00          -1.272620D-01
      5.215000D-01           5.445290D-01
      1.596000D-01           5.804960D-01
S   1   1.00
      1.596000D-01           1.000000D+00
P   4   1.00
      9.439000D+00           3.810900D-02
      2.002000D+00           2.094800D-01
      5.456000D-01           5.085570D-01
      1.517000D-01           4.688420D-01
P   1   1.00
      1.517000D-01           1.000000D+00
D   1   1.00
      5.500000D-01           1.0000000
****
N     0
S   9   1.00
      9.046000D+03           7.000000D-04
      1.357000D+03           5.389000D-03
      3.093000D+02           2.740600D-02
      8.773000D+01           1.032070D-01
      2.856000D+01           2.787230D-01
      1.021000D+01           4.485400D-01
      3.838000D+00           2.782380D-01
      7.466000D-01           1.544000D-02
      2.248000D-01          -2.864000D-03
S   9   1.00
      9.046000D+03          -1.530000D-04
      1.357000D+03          -1.208000D-03
      3.093000D+02          -5.992000D-03
      8.773000D+01          -2.454400D-02
      2.856000D+01          -6.745900D-02
      1.021000D+01          -1.580780D-01
      3.838000D+00          -1.218310D-01
      7.466000D-01           5.490030D-01
      2.248000D-01           5.788150D-01
S   1   1.00
      2.248000D-01           1.000000D+00
P   4   1.00
      1.355000D+01           3.991900D-02
      2.917000D+00           2.171690D-01
      7.973000D-01           5.103190D-01
      2.185000D-01           4.622140D-01
P   1   1.00
      2.185000D-01           1.000000D+00
D   1   1.00
      8.170000D-01           1.0000000
****
O     0
S   9   1.00
      1.172000D+04           7.100000D-04
      1.759000D+03           5.470000D-03
      4.008000D+02           2.783700D-02
      1.137000D+02           1.048000D-01
      3.703000D+01           2.830620D-01
      1.327000D+01           4.487190D-01
      5.025000D+00           2.709520D-01
      1.013000D+00           1.545800D-02
      3.023000D-01          -2.585000D-03
S   9   1.00
      1.172000D+04          -1.600000D-04
      1.759000D+03          -1.263000D-03
      4.008000D+02          -6.267000D-03
      1.137000D+02          -2.571600D-02
      3.703000D+01          -7.092400D-02
      1.327000D+01          -1.654110D-01
      5.025000D+00          -1.169550D-01
      1.013000D+00           5.573680D-01
      3.023000D-01           5.727590D-01
S   1   1.00
      3.023000D-01           1.000000D+00
P   4   1.00
      1.770000D+01           4.301800D-02
      3.854000D+00           2.289130D-01
      1.046000D+00           5.087280D-01
      2.753000D-01           4.605310D-01
P   1   1.00
      2.753000D-01           1.000000D+00
D   1   1.00
      1.185000D+00           1.0000000
****
S     0
S   12   1.00
      1.108000D+05           2.476350D-04
      1.661000D+04           1.920260D-03
      3.781000D+03           9.961920D-03
      1.071000D+03           4.029750D-02
      3.498000D+02           1.286040D-01
      1.263000D+02           3.034800D-01
      4.926000D+01           4.214320D-01
      2.016000D+01           2.307810D-01
      5.720000D+00           1.789710D-02
      2.182000D+00          -2.975160D-03
      4.327000D-01           8.495220D-04
      1.570000D-01          -3.679360D-04
S   12   1.00
      1.108000D+05          -6.870390D-05
      1.661000D+04          -5.276810D-04
      3.781000D+03          -2.796710D-03
      1.071000D+03          -1.126510D-02
      3.498000D+02          -3.888340D-02
      1.263000D+02          -9.950250D-02
      4.926000D+01          -1.997400D-01
      2.016000D+01          -1.233600D-01
      5.720000D+00           5.131940D-01
      2.182000D+00           6.071200D-01
      4.327000D-01           3.967530D-02
      1.570000D-01          -9.468640D-03
S   12   1.00
      1.108000D+05           1.990770D-05
      1.661000D+04           1.534830D-04
      3.781000D+03           8.095030D-04
      1.071000D+03           3.289740D-03
      3.498000D+02           1.129670D-02
      1.263000D+02           2.963850D-02
      4.926000D+01           5.998510D-02
      2.016000D+01           4.132480D-02
      5.720000D+00          -2.074740D-01
      2.182000D+00          -3.928890D-01
      4.327000D-01           6.328400D-01
      1.570000D-01           5.569240D-01
S   1   1.00
      1.570000D-01           1.000000D+00
P   8   1.00
      3.997000D+02           4.475410D-03
      9.419000D+01           3.417080D-02
      2.975000D+01           1.442500D-01
      1.077000D+01           3.539280D-01
      4.119000D+00           4.590850D-01
      1.625000D+00           2.063830D-01
      4.726000D-01           1.021410D-02
      1.407000D-01          -6.031220D-05
P   8   1.00
      3.997000D+02          -1.162510D-03
      9.419000D+01          -8.656640D-03
      2.975000D+01          -3.908860D-02
      1.077000D+01          -9.346250D-02
      4.119000D+00          -1.479940D-01
      1.625000D+00           3.019040D-02
      4.726000D-01           5.615730D-01
      1.407000D-01           5.347760D-01
P   1   1.00
      1.407000D-01           1.000000D+00
D   1   1.00
      4.790000D-01           1.0000000
****
*/
BasisSet<double> make_ccpvdz( const Molecule& mol, SphericalType sph ) {

#if 1
  std::string basis_path = GAUXC_REF_DATA_PATH  "/../basis/old/cc-pvdz.g94";
  return parse_basis( mol, basis_path, sph );
#else
  BasisSet<double> basis;

  for( const auto& atom : mol ) {


    if( atom.Z == 1ll ) {

      basis.emplace_back(
        Shell<double>(
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          { 1.301000e+01, 1.962000e+00, 4.446000e-01,
            1.220000e-01 },
          { 1.968500e-02, 1.379770e-01, 4.781480e-01,
            5.012400e-01 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {1.220000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          {7.270000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z == 6ll ) {

      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 6.665000e+03, 1.000000e+03, 2.280000e+02,
            6.471000e+01, 2.106000e+01, 7.495000e+00,
            2.797000e+00, 5.215000e-01, 1.596000e-01 },
          { 6.920000e-04, 5.329000e-03,  2.707700e-02,
            1.017180e-01, 2.747400e-01,  4.485640e-01,
            2.850740e-01, 1.520400e-02, -3.191000e-03 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 6.665000e+03, 1.000000e+03, 2.280000e+02,
            6.471000e+01, 2.106000e+01, 7.495000e+00,
            2.797000e+00, 5.215000e-01, 1.596000e-01 },
          { -1.460000e-04, -1.154000e-03, -5.725000e-03,
            -2.331200e-02, -6.395500e-02, -1.499810e-01,
            -1.272620e-01,  5.445290e-01,  5.804960e-01 },
          { atom.x, atom.y, atom.z }
        )
      );


      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {1.596000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 9.439000e+00, 2.002000e+00, 5.456000e-01,
            1.517000e-01 }, 
          { 3.810900e-02, 2.094800e-01, 5.085570e-01,
            4.688420e-01 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          {1.517000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(2), SphericalType(sph),
          {5.500000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z == 7ll ) {

      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 9.046000e+03, 1.357000e+03, 3.093000e+02,
            8.773000e+01, 2.856000e+01, 1.021000e+01,
            3.838000e+00, 7.466000e-01, 2.248000e-01 },
          { 7.000000e-04, 5.389000e-03, 2.740600e-02,
            1.032070e-01, 2.787230e-01, 4.485400e-01,
            2.782380e-01, 1.544000e-02, -2.864000e-03},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 9.046000e+03, 1.357000e+03, 3.093000e+02,
            8.773000e+01, 2.856000e+01, 1.021000e+01,
            3.838000e+00, 7.466000e-01, 2.248000e-01 },
          { -1.530000e-04, -1.208000e-03, -5.992000e-03,
            -2.454400e-02, -6.745900e-02, -1.580780e-01,
            -1.218310e-01,  5.490030e-01,  5.788150e-01 },
          { atom.x, atom.y, atom.z }
        )
      );
      
      
      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {2.248000e-01}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 1.355000e+01, 2.917000e+00, 7.973000e-01,
            2.185000e-01 },
          { 3.991900e-02, 2.171690e-01, 5.103190e-01,
            4.622140e-01 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          {2.185000e-01}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(2), SphericalType(sph),
          {8.170000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );


    } else if( atom.Z == 8ll ) {




      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 1.172000e+04, 1.759000e+03, 4.008000e+02,
            1.137000e+02, 3.703000e+01, 1.327000e+01,
            5.025000e+00, 1.013000e+00, 3.023000e-01 },
          { 7.100000e-04, 5.470000e-03,  2.783700e-02,
            1.048000e-01, 2.830620e-01,  4.487190e-01,
            2.709520e-01, 1.545800e-02, -2.585000e-03 }, 
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(8), AngularMomentum(0), SphericalType(false),
          { 1.172000e+04, 1.759000e+03, 4.008000e+02,
            1.137000e+02, 3.703000e+01, 1.327000e+01,
            5.025000e+00, 1.013000e+00, 3.023000e-01 },
          { -1.600000e-04, -1.263000e-03, -6.267000e-03,
            -2.571600e-02, -7.092400e-02, -1.654110e-01,
            -1.169550e-01,  5.573680e-01,  5.727590e-01 },
          { atom.x, atom.y, atom.z }
        )
      );


      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {3.023000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );



      basis.emplace_back(
        Shell<double>(
          PrimSize(3), AngularMomentum(1), SphericalType(false),
          { 1.770000e+01, 3.854000e+00, 1.046000e+00,
            2.753000e-01 },
          { 4.301800e-02, 2.289130e-01, 5.087280e-01,
            4.605310e-01 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          {2.753000e-1}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(2), SphericalType(sph),
          {1.185000e+0}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else if( atom.Z == 16ll ) {

      basis.emplace_back(
        Shell<double>(
          PrimSize(11), AngularMomentum(0), SphericalType(false),
          { 1.108000e+05, 1.661000e+04, 3.781000e+03,
            1.071000e+03, 3.498000e+02, 1.263000e+02,
            4.926000e+01, 2.016000e+01, 5.720000e+00,
            2.182000e+00, 4.327000e-01, 1.570000e-01 }, 
          { 2.476350e-04, 1.920260e-03, 9.961920e-03,
            4.029750e-02, 1.286040e-01, 3.034800e-01,
            4.214320e-01, 2.307810e-01, 1.789710e-02,
           -2.975160e-03, 8.495220e-04, -3.679360e-04}, 
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(11), AngularMomentum(0), SphericalType(false),
          { 1.108000e+05, 1.661000e+04, 3.781000e+03,
            1.071000e+03, 3.498000e+02, 1.263000e+02,
            4.926000e+01, 2.016000e+01, 5.720000e+00,
            2.182000e+00, 4.327000e-01, 1.570000e-01 }, 
          {-6.870390e-05, -5.276810e-04, -2.796710e-03,
           -1.126510e-02, -3.888340e-02, -9.950250e-02,
           -1.997400e-01, -1.233600e-01, 5.131940e-01,
            6.071200e-01, 3.967530e-02, -9.468640e-03  }, 
          { atom.x, atom.y, atom.z }
        )
      );
      
      basis.emplace_back(
        Shell<double>(
          PrimSize(11), AngularMomentum(0), SphericalType(false),
          { 1.108000e+05, 1.661000e+04, 3.781000e+03,
            1.071000e+03, 3.498000e+02, 1.263000e+02,
            4.926000e+01, 2.016000e+01, 5.720000e+00,
            2.182000e+00, 4.327000e-01, 1.570000e-01 }, 
          { 1.990770e-05, 1.534830e-04, 8.095030e-04,
            3.289740e-03, 1.129670e-02, 2.963850e-02,
            5.998510e-02, 4.132480e-02, -2.074740e-01,
           -3.928890e-01, 6.328400e-01, 5.569240e-01  }, 
          { atom.x, atom.y, atom.z }
        )
      );
      
      
      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(0), SphericalType(false),
          {1.570000e-01}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );
      
            
      basis.emplace_back(
        Shell<double>(
          PrimSize(7), AngularMomentum(1), SphericalType(false),
          { 3.997000e+02, 9.419000e+01, 2.975000e+01,
            1.077000e+01, 4.119000e+00, 1.625000e+00,
            4.726000e-01, 1.407000e-01 }, 
          { 4.475410e-03, 3.417080e-02, 1.442500e-01,
            3.539280e-01, 4.590850e-01, 2.063830e-01,
            1.021410e-02, -6.031220e-05 },
          { atom.x, atom.y, atom.z }
        )
      );

      basis.emplace_back(
        Shell<double>(
          PrimSize(7), AngularMomentum(1), SphericalType(false),
          { 3.997000e+02, 9.419000e+01, 2.975000e+01,
            1.077000e+01, 4.119000e+00, 1.625000e+00,
            4.726000e-01, 1.407000e-01 }, 
          { -1.162510e-03, -8.656640e-03, -3.908860e-02,
            -9.346250e-02, -1.479940e-01, 3.019040e-02,
             5.615730e-01, 5.347760e-01 },
          { atom.x, atom.y, atom.z }
        )
      );
       
       
      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(1), SphericalType(false),
          {1.407000e-01}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );
       
      basis.emplace_back(
        Shell<double>(
          PrimSize(1), AngularMomentum(2), SphericalType(sph),
          {4.790000e-01}, {1.000000},
          { atom.x, atom.y, atom.z }
        )
      );

    } else { throw std::runtime_error("NYI"); }


  }

  return basis;
#endif
}




}
