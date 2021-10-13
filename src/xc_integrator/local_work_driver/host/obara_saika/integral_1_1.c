void integral_1_1(int npts,
                  shell shellA,
                  shell shellB,
                  point *points,
                  double *Xi,
                  int ldX,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   double temp[9];

   for(int point_idx = 0; point_idx < nr_points; ++point_idx) {
      point C = *(point_list + point_idx);

      double xA = shellA.origin.x;
      double yA = shellA.origin.y;
      double zA = shellA.origin.z;

      double xB = shellB.origin.x;
      double yB = shellB.origin.y;
      double zB = shellB.origin.z;

      double X_AB = (xA - xB);
      double Y_AB = (yA - yB);
      double Z_AB = (zA - zB);

      double beta_in = 0.0;
      for(int i = 0; i < shellA.m; ++i) {
         for(int j = 0; j < shellB.m; ++j) {
            double aA = shellA.coeff[i].alpha;
            double cA = shellA.coeff[i].coeff;

            double aB = shellB.coeff[j].alpha;
            double cB = shellB.coeff[j].coeff;

            double RHO = aA + aB;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double xP = (aA * xA + aB * xB) * RHO_INV;
            double yP = (aA * yA + aB * yB) * RHO_INV;
            double zP = (aA * zA + aB * zB) * RHO_INV;

            double X_PB = (xP - xB);
            double Y_PB = (yP - yB);
            double Z_PB = (zP - zB);

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t10, t11, t20;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);

            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            *(temp + 0) = beta_in * (*(temp + 0)) + t10;

            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 3) = beta_in * (*(temp + 3)) + t20;

            t20 = Y_PA * t10 - Y_PC * t11;
            *(temp + 4) = beta_in * (*(temp + 4)) + t20;

            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 5) = beta_in * (*(temp + 5)) + t20;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            *(temp + 1) = beta_in * (*(temp + 1)) + t10;

            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 6) = beta_in * (*(temp + 6)) + t20;

            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 7) = beta_in * (*(temp + 7)) + t20;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            *(temp + 2) = beta_in * (*(temp + 2)) + t10;

            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 8) = beta_in * (*(temp + 8)) + t20;

            beta_in = 1.0;
         }
      }

      double *Xik = *(Xi + point_idx * ldX);
      double *Gjk = *(Gj + point_idx * ldG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 3) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      t1 = *(temp + 4) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      t2 = *(temp + 5) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;

      X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 4) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 6) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 7) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;

      Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 5) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 7) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 8) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;

      Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;

      *(Gjk + 0) *= *(weights + point_idx);
      *(Gjk + 1) *= *(weights + point_idx);
      *(Gjk + 2) *= *(weights + point_idx);
   }
}
