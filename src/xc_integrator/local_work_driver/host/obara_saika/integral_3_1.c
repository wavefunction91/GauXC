void integral_3_1(int npts,
                  shell shellA,
                  shell shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   double temp[25];

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

            double t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);
            t03 = eval * boys_reference(3, tval);
            t04 = eval * boys_reference(4, tval);

            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t13 = X_PA * t03 - X_PC * t04;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PA * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 0) = beta_in * (*(temp + 0)) + t30;

            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 10) = beta_in * (*(temp + 10)) + t40;

            t40 = Y_PA * t30 - Y_PC * t31;
            *(temp + 11) = beta_in * (*(temp + 11)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 12) = beta_in * (*(temp + 12)) + t40;

            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            *(temp + 1) = beta_in * (*(temp + 1)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 13) = beta_in * (*(temp + 13)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 14) = beta_in * (*(temp + 14)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 2) = beta_in * (*(temp + 2)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 15) = beta_in * (*(temp + 15)) + t40;

            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 3) = beta_in * (*(temp + 3)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 16) = beta_in * (*(temp + 16)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 17) = beta_in * (*(temp + 17)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 4) = beta_in * (*(temp + 4)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 18) = beta_in * (*(temp + 18)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 5) = beta_in * (*(temp + 5)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 19) = beta_in * (*(temp + 19)) + t40;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 6) = beta_in * (*(temp + 6)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 20) = beta_in * (*(temp + 20)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 21) = beta_in * (*(temp + 21)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 7) = beta_in * (*(temp + 7)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 22) = beta_in * (*(temp + 22)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 8) = beta_in * (*(temp + 8)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 23) = beta_in * (*(temp + 23)) + t40;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 9) = beta_in * (*(temp + 9)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 24) = beta_in * (*(temp + 24)) + t40;

            beta_in = 1.0;
         }
      }

      double *Xik = *(Xi + point_idx * ldX);
      double *Xjk = *(Xj + point_idx * ldX);
      double *Gik = *(Gi + point_idx * ldG);
      double *Gjk = *(Gj + point_idx * ldG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 10) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      *(Gik + 0) += *(Xjk + 0) * t0;
      t1 = *(temp + 11) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      *(Gik + 0) += *(Xjk + 1) * t1;
      t2 = *(temp + 12) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      *(Gik + 0) += *(Xjk + 2) * t2;
      t3 = *(temp + 13) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      *(Gik + 0) += *(Xjk + 3) * t3;
      t4 = *(temp + 14) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      *(Gik + 0) += *(Xjk + 4) * t4;
      t5 = *(temp + 15) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;
      *(Gik + 0) += *(Xjk + 5) * t5;
      t6 = *(temp + 16) * const_value;
      *(Gjk + 6) += *(Xik + 0) * t6;
      *(Gik + 0) += *(Xjk + 6) * t6;
      t7 = *(temp + 17) * const_value;
      *(Gjk + 7) += *(Xik + 0) * t7;
      *(Gik + 0) += *(Xjk + 7) * t7;
      t8 = *(temp + 18) * const_value;
      *(Gjk + 8) += *(Xik + 0) * t8;
      *(Gik + 0) += *(Xjk + 8) * t8;
      t9 = *(temp + 19) * const_value;
      *(Gjk + 9) += *(Xik + 0) * t9;
      *(Gik + 0) += *(Xjk + 9) * t9;

      X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      *(Gik + 0) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      *(Gik + 0) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      *(Gik + 0) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      *(Gik + 0) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      *(Gik + 0) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;
      *(Gik + 0) += *(Xjk + 5) * t5;
      t6 = *(temp + 6) * const_value;
      *(Gjk + 6) += *(Xik + 0) * t6;
      *(Gik + 0) += *(Xjk + 6) * t6;
      t7 = *(temp + 7) * const_value;
      *(Gjk + 7) += *(Xik + 0) * t7;
      *(Gik + 0) += *(Xjk + 7) * t7;
      t8 = *(temp + 8) * const_value;
      *(Gjk + 8) += *(Xik + 0) * t8;
      *(Gik + 0) += *(Xjk + 8) * t8;
      t9 = *(temp + 9) * const_value;
      *(Gjk + 9) += *(Xik + 0) * t9;
      *(Gik + 0) += *(Xjk + 9) * t9;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 11) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      *(Gik + 1) += *(Xjk + 0) * t0;
      t1 = *(temp + 13) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      *(Gik + 1) += *(Xjk + 1) * t1;
      t2 = *(temp + 14) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      *(Gik + 1) += *(Xjk + 2) * t2;
      t3 = *(temp + 16) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      *(Gik + 1) += *(Xjk + 3) * t3;
      t4 = *(temp + 17) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      *(Gik + 1) += *(Xjk + 4) * t4;
      t5 = *(temp + 18) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;
      *(Gik + 1) += *(Xjk + 5) * t5;
      t6 = *(temp + 20) * const_value;
      *(Gjk + 6) += *(Xik + 1) * t6;
      *(Gik + 1) += *(Xjk + 6) * t6;
      t7 = *(temp + 21) * const_value;
      *(Gjk + 7) += *(Xik + 1) * t7;
      *(Gik + 1) += *(Xjk + 7) * t7;
      t8 = *(temp + 22) * const_value;
      *(Gjk + 8) += *(Xik + 1) * t8;
      *(Gik + 1) += *(Xjk + 8) * t8;
      t9 = *(temp + 23) * const_value;
      *(Gjk + 9) += *(Xik + 1) * t9;
      *(Gik + 1) += *(Xjk + 9) * t9;

      Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      *(Gik + 1) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      *(Gik + 1) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      *(Gik + 1) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      *(Gik + 1) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      *(Gik + 1) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;
      *(Gik + 1) += *(Xjk + 5) * t5;
      t6 = *(temp + 6) * const_value;
      *(Gjk + 6) += *(Xik + 1) * t6;
      *(Gik + 1) += *(Xjk + 6) * t6;
      t7 = *(temp + 7) * const_value;
      *(Gjk + 7) += *(Xik + 1) * t7;
      *(Gik + 1) += *(Xjk + 7) * t7;
      t8 = *(temp + 8) * const_value;
      *(Gjk + 8) += *(Xik + 1) * t8;
      *(Gik + 1) += *(Xjk + 8) * t8;
      t9 = *(temp + 9) * const_value;
      *(Gjk + 9) += *(Xik + 1) * t9;
      *(Gik + 1) += *(Xjk + 9) * t9;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 12) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      *(Gik + 2) += *(Xjk + 0) * t0;
      t1 = *(temp + 14) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      *(Gik + 2) += *(Xjk + 1) * t1;
      t2 = *(temp + 15) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      *(Gik + 2) += *(Xjk + 2) * t2;
      t3 = *(temp + 17) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      *(Gik + 2) += *(Xjk + 3) * t3;
      t4 = *(temp + 18) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      *(Gik + 2) += *(Xjk + 4) * t4;
      t5 = *(temp + 19) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;
      *(Gik + 2) += *(Xjk + 5) * t5;
      t6 = *(temp + 21) * const_value;
      *(Gjk + 6) += *(Xik + 2) * t6;
      *(Gik + 2) += *(Xjk + 6) * t6;
      t7 = *(temp + 22) * const_value;
      *(Gjk + 7) += *(Xik + 2) * t7;
      *(Gik + 2) += *(Xjk + 7) * t7;
      t8 = *(temp + 23) * const_value;
      *(Gjk + 8) += *(Xik + 2) * t8;
      *(Gik + 2) += *(Xjk + 8) * t8;
      t9 = *(temp + 24) * const_value;
      *(Gjk + 9) += *(Xik + 2) * t9;
      *(Gik + 2) += *(Xjk + 9) * t9;

      Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      *(Gik + 2) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      *(Gik + 2) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      *(Gik + 2) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      *(Gik + 2) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      *(Gik + 2) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;
      *(Gik + 2) += *(Xjk + 5) * t5;
      t6 = *(temp + 6) * const_value;
      *(Gjk + 6) += *(Xik + 2) * t6;
      *(Gik + 2) += *(Xjk + 6) * t6;
      t7 = *(temp + 7) * const_value;
      *(Gjk + 7) += *(Xik + 2) * t7;
      *(Gik + 2) += *(Xjk + 7) * t7;
      t8 = *(temp + 8) * const_value;
      *(Gjk + 8) += *(Xik + 2) * t8;
      *(Gik + 2) += *(Xjk + 8) * t8;
      t9 = *(temp + 9) * const_value;
      *(Gjk + 9) += *(Xik + 2) * t9;
      *(Gik + 2) += *(Xjk + 9) * t9;

      *(Gjk + 0) *= *(weights + point_idx);
      *(Gjk + 1) *= *(weights + point_idx);
      *(Gjk + 2) *= *(weights + point_idx);
      *(Gjk + 3) *= *(weights + point_idx);
      *(Gjk + 4) *= *(weights + point_idx);
      *(Gjk + 5) *= *(weights + point_idx);
      *(Gjk + 6) *= *(weights + point_idx);
      *(Gjk + 7) *= *(weights + point_idx);
      *(Gjk + 8) *= *(weights + point_idx);
      *(Gjk + 9) *= *(weights + point_idx);

      *(Gik + 0) *= *(weights + point_idx);
      *(Gik + 1) *= *(weights + point_idx);
      *(Gik + 2) *= *(weights + point_idx);
   }
}
