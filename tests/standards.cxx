#include "standards.hpp"

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

  return mol;
}

Molecule make_ubiquitin() {

  Molecule mol;

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

  return mol;
}



/*
H     0
S   3   1.00
      0.1873113696D+02       0.3349460434D-01
      0.2825394365D+01       0.2347269535D+00
      0.6401216923D+00       0.8137573261D+00
S   1   1.00
      0.1612777588D+00       1.0000000
****
C     0
S   6   1.00
      0.3047524880D+04       0.1834737132D-02
      0.4573695180D+03       0.1403732281D-01
      0.1039486850D+03       0.6884262226D-01
      0.2921015530D+02       0.2321844432D+00
      0.9286662960D+01       0.4679413484D+00
      0.3163926960D+01       0.3623119853D+00
SP   3   1.00
      0.7868272350D+01      -0.1193324198D+00       0.6899906659D-01
      0.1881288540D+01      -0.1608541517D+00       0.3164239610D+00
      0.5442492580D+00       0.1143456438D+01       0.7443082909D+00
SP   1   1.00
      0.1687144782D+00       0.1000000000D+01       0.1000000000D+01
D   1   1.00
      0.8000000000D+00       1.0000000
****
N     0
S   6   1.00
      0.4173511460D+04       0.1834772160D-02
      0.6274579110D+03       0.1399462700D-01
      0.1429020930D+03       0.6858655181D-01
      0.4023432930D+02       0.2322408730D+00
      0.1282021290D+02       0.4690699481D+00
      0.4390437010D+01       0.3604551991D+00
SP   3   1.00
      0.1162636186D+02      -0.1149611817D+00       0.6757974388D-01
      0.2716279807D+01      -0.1691174786D+00       0.3239072959D+00
      0.7722183966D+00       0.1145851947D+01       0.7408951398D+00
SP   1   1.00
      0.2120314975D+00       0.1000000000D+01       0.1000000000D+01
D   1   1.00
      0.8000000000D+00       1.0000000
****
O     0
S   6   1.00
      0.5484671660D+04       0.1831074430D-02
      0.8252349460D+03       0.1395017220D-01
      0.1880469580D+03       0.6844507810D-01
      0.5296450000D+02       0.2327143360D+00
      0.1689757040D+02       0.4701928980D+00
      0.5799635340D+01       0.3585208530D+00
SP   3   1.00
      0.1553961625D+02      -0.1107775495D+00       0.7087426823D-01
      0.3599933586D+01      -0.1480262627D+00       0.3397528391D+00
      0.1013761750D+01       0.1130767015D+01       0.7271585773D+00
SP   1   1.00
      0.2700058226D+00       0.1000000000D+01       0.1000000000D+01
D   1   1.00
      0.8000000000D+00       1.0000000
****
S     0
S   6   1.00
      0.2191710000D+05       0.1869240849D-02
      0.3301490000D+04       0.1423030646D-01
      0.7541460000D+03       0.6969623166D-01
      0.2127110000D+03       0.2384871083D+00
      0.6798960000D+02       0.4833072195D+00
      0.2305150000D+02       0.3380741536D+00
SP   6   1.00
      0.4237350000D+03      -0.2376770499D-02       0.4061009982D-02
      0.1007100000D+03      -0.3169300665D-01       0.3068129986D-01
      0.3215990000D+02      -0.1133170238D+00       0.1304519994D+00
      0.1180790000D+02       0.5609001177D-01       0.3272049985D+00
      0.4631100000D+01       0.5922551243D+00       0.4528509980D+00
      0.1870250000D+01       0.4550060955D+00       0.2560419989D+00
SP   3   1.00
      0.2615840000D+01      -0.2503731142D+00      -0.1451048955D-01
      0.9221670000D+00       0.6695676310D-01       0.3102627765D+00
      0.3412870000D+00       0.1054506269D+01       0.7544824565D+00
SP   1   1.00
      0.1171670000D+00       0.1000000000D+01       0.1000000000D+01
D   1   1.00
      0.6500000000D+00       1.0000000
****
*/
BasisSet<double> make_631Gd( const Molecule& mol, SphericalType sph ) {

  BasisSet<double> basis;

  for( const auto& atom : mol ) {

    if( atom.Z == 1ll ) {

      basis.emplace_back(
        Shell<double>( 
          PrimSize(3), AngularMomentum(0), SphericalType(false),
          {0.3349460434e-01, 0.2347269535e+00, 0.8137573261e+00},
          {0.1873113696e+02, 0.2825394365e+01, 0.6401216923e+00},
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
          {1.185000e-1}, {1.000000},
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
}




}
