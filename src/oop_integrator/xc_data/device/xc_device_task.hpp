#pragma once

namespace GauXC {

struct XCDeviceTask {

  size_t nbe;
  size_t npts;
  size_t ncut;
  size_t nblock;
  size_t nshells;

  double* points;
  double* weights;
  size_t* shell_list;
  size_t* shell_offs;
  int32_t* submat_cut;
  int32_t* submat_block;

  double*   nbe_scr;
  double*   zmat;
  double*   bf, *dbfx, *dbfy, *dbfz;
  double*   den, *ddenx, *ddeny, *ddenz;
  double*   eps, *gamma;
  double*   vrho, *vgamma;

  size_t iParent;
  double dist_nearest;
  double * dist_scratch;

};

}
