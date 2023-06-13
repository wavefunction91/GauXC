
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct node {
  int iA, jA, kA;
  int iB, jB, kB;
  
  int level;
  int vars;

  int valid;
  int offset;

  char var_pa[5];
  char var_pc[5];
  
  int nr_children;
  struct node *children[3];
};

void traverseX_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list);
void traverseY_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list);
void traverseZ_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list);

void traverseX_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list) {
  int level = iA + jA + kA;
  int offset = offset_list[level];
  int vars = lA + lB - level;
  
  node_list[offset].iA = iA;
  node_list[offset].jA = jA;
  node_list[offset].kA = kA;

  node_list[offset].iB = 0;
  node_list[offset].jB = 0;
  node_list[offset].kB = 0;

  node_list[offset].level = level;
  node_list[offset].vars = vars + 1;

  node_list[offset].valid = (level >= lA) ? 1 : 0;
  node_list[offset].offset = offset - partial_size;

  node_list[offset].nr_children = 0;
  node_list[offset].children[0] = NULL;
  node_list[offset].children[1] = NULL;
  node_list[offset].children[2] = NULL;

  sprintf(node_list[offset].var_pa, "X_PA");
  sprintf(node_list[offset].var_pc, "X_PC");
  
  offset_list[level]++;
  
  if(vars > 0) {
    node_list[offset].nr_children = 3;
    
    node_list[offset].children[0] = &(node_list[offset_list[level + 1]]);
    traverseX_init_dfs(iA + 1, jA, kA, lA, lB, partial_size, node_list, offset_list);

    node_list[offset].children[1] = &(node_list[offset_list[level + 1]]);
    traverseY_init_dfs(iA, jA + 1, kA, lA, lB, partial_size, node_list, offset_list);

    node_list[offset].children[2] = &(node_list[offset_list[level + 1]]); 
    traverseZ_init_dfs(iA, jA, kA + 1, lA, lB, partial_size, node_list, offset_list);
  }
}

void traverseY_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list) {
  int level = iA + jA + kA;
  int offset = offset_list[level];
  int vars = lA + lB - level;
  
  node_list[offset].iA = iA;
  node_list[offset].jA = jA;
  node_list[offset].kA = kA;

  node_list[offset].iB = 0;
  node_list[offset].jB = 0;
  node_list[offset].kB = 0;

  node_list[offset].level = level;
  node_list[offset].vars = vars + 1;

  node_list[offset].valid = (level >= lA) ? 1 : 0;
  node_list[offset].offset = offset - partial_size;

  node_list[offset].nr_children = 0;
  node_list[offset].children[0] = NULL;
  node_list[offset].children[1] = NULL;
  node_list[offset].children[2] = NULL;

  sprintf(node_list[offset].var_pa, "Y_PA");
  sprintf(node_list[offset].var_pc, "Y_PC");
  
  offset_list[level]++;
  
  if(vars > 0) {
    node_list[offset].nr_children = 2;
    
    node_list[offset].children[0] = &(node_list[offset_list[level + 1]]);
    traverseY_init_dfs(iA, jA + 1, kA, lA, lB, partial_size, node_list, offset_list);

    node_list[offset].children[1] = &(node_list[offset_list[level + 1]]); 
    traverseZ_init_dfs(iA, jA, kA + 1, lA, lB, partial_size, node_list, offset_list);
  }
}

void traverseZ_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list) {
  int level = iA + jA + kA;
  int offset = offset_list[level];
  int vars = lA + lB - level;
  
  node_list[offset].iA = iA;
  node_list[offset].jA = jA;
  node_list[offset].kA = kA;

  node_list[offset].iB = 0;
  node_list[offset].jB = 0;
  node_list[offset].kB = 0;

  node_list[offset].level = level;
  node_list[offset].vars = vars + 1;

  node_list[offset].valid = (level >= lA) ? 1 : 0;
  node_list[offset].offset = offset - partial_size;

  node_list[offset].nr_children = 0;
  node_list[offset].children[0] = NULL;
  node_list[offset].children[1] = NULL;
  node_list[offset].children[2] = NULL;

  sprintf(node_list[offset].var_pa, "Z_PA");
  sprintf(node_list[offset].var_pc, "Z_PC");
  
  offset_list[level]++;
  
  if(vars > 0) {
    node_list[offset].nr_children = 1;
 
    node_list[offset].children[0] = &(node_list[offset_list[level + 1]]); 
    traverseZ_init_dfs(iA, jA, kA + 1, lA, lB, partial_size, node_list, offset_list);
  }
}

void traverse_init_dfs(int iA, int jA, int kA, int lA, int lB, int partial_size, struct node *node_list, int *offset_list) {
  int level = iA + jA + kA;
  int offset = offset_list[level];
  int vars = lA + lB - level;
  
  node_list[offset].iA = iA;
  node_list[offset].jA = jA;
  node_list[offset].kA = kA;

  node_list[offset].iB = 0;
  node_list[offset].jB = 0;
  node_list[offset].kB = 0;

  node_list[offset].level = level;
  node_list[offset].vars = vars + 1;

  node_list[offset].valid = (level >= lA) ? 1 : 0;
  node_list[offset].offset = offset - partial_size;

  node_list[offset].nr_children = 0;
  node_list[offset].children[0] = NULL;
  node_list[offset].children[1] = NULL;
  node_list[offset].children[2] = NULL;
  
  offset_list[level]++;
  
  if(vars > 0) {
    node_list[offset].nr_children = 3;
    
    node_list[offset].children[0] = &(node_list[offset_list[level + 1]]);
    traverseX_init_dfs(iA + 1, jA, kA, lA, lB, partial_size, node_list, offset_list);

    node_list[offset].children[1] = &(node_list[offset_list[level + 1]]);
    traverseY_init_dfs(iA, jA + 1, kA, lA, lB, partial_size, node_list, offset_list);

    node_list[offset].children[2] = &(node_list[offset_list[level + 1]]); 
    traverseZ_init_dfs(iA, jA, kA + 1, lA, lB, partial_size, node_list, offset_list);
  }
}

void initialize_tree_structure(int type, int lA, int lB, int size, struct node *node_list) {
  int partial_size = 0;
  for(int i = 0; i < lA; ++i) {
    partial_size += (i + 1) * (i + 2) / 2;
  }
  
  int *offset_list = (int*) malloc((lA + lB + 1) * sizeof(int));

  int offset = 0;
  for(int i = 0; i < lA + lB + 1; ++i) {
    offset_list[i] = offset;
    offset += (i + 1) * (i + 2) / 2;
  }

  // initialization part
  traverse_init_dfs(0, 0, 0, lA, lB, partial_size, node_list, offset_list);
  
  free(offset_list);
}

void traverse_dfs_vrr(FILE *f, int lA, int lB, struct node *root_node) {
  if(root_node != NULL) {
    if(root_node -> level == 0) {
      for(int v = 0; v < root_node -> vars; ++v) {
	fprintf(f, "         t%d%d = SCALAR_MUL(eval, t%d%d);\n", root_node -> level, v, root_node -> level, v);
      }
    } else if (root_node -> level == 1) {
      for(int v = 0; v < root_node -> vars; ++v) {
	fprintf(f, "         t%d%d = SCALAR_MUL(%s, t%d%d);\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v);
	fprintf(f, "         t%d%d = SCALAR_FNMA(%s, t%d%d, t%d%d);\n", root_node -> level, v, root_node -> var_pc, root_node -> level - 1, v + 1, root_node -> level, v);
      }
    } else {
      int iteration = 0;
      if(strcmp(root_node -> var_pa, "X_PA") == 0) {
	iteration = root_node -> iA - 1;
      } else if(strcmp(root_node -> var_pa, "Y_PA") == 0) {
	iteration = root_node -> jA - 1;
      } else {
	iteration = root_node -> kA - 1;
      }

      if(iteration == 0) {
	for(int v = 0; v < root_node -> vars; ++v) {
	  fprintf(f, "         t%d%d = SCALAR_MUL(%s, t%d%d);\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v);
	  fprintf(f, "         t%d%d = SCALAR_FNMA(%s, t%d%d, t%d%d);\n", root_node -> level, v, root_node -> var_pc, root_node -> level - 1, v + 1, root_node -> level, v);
	}
      } else {
	for(int v = 0; v < root_node -> vars; ++v) {
	  fprintf(f, "         t%d%d = SCALAR_MUL(%s, t%d%d);\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v);
	  fprintf(f, "         t%d%d = SCALAR_FNMA(%s, t%d%d, t%d%d);\n", root_node -> level, v, root_node -> var_pc, root_node -> level - 1, v + 1, root_node -> level, v);
	  fprintf(f, "         tx = SCALAR_SUB(t%d%d, t%d%d);\n", root_node -> level - 2, v, root_node ->level - 2, v + 1);
	  fprintf(f, "         ty = SCALAR_SET1(0.5 * %d);\n", iteration);
	  fprintf(f, "         ty = SCALAR_MUL(ty, RHO_INV);\n");
	  fprintf(f, "         t%d%d = SCALAR_FMA(tx, ty, t%d%d);\n", root_node -> level, v, root_node -> level, v);
	}
      }
    }

    if(root_node -> valid) {
      fprintf(f, "         tx = SCALAR_LOAD((temp + %d * blockDim.x + threadIdx.x));\n", root_node -> offset);
      fprintf(f, "         tx = SCALAR_ADD(tx, t%d%d);\n", root_node -> level, 0);
      fprintf(f, "         SCALAR_STORE((temp + %d * blockDim.x + threadIdx.x), tx);\n", root_node -> offset);
    }
    
    for(int i = 0; i < root_node -> nr_children; ++i) {
      traverse_dfs_vrr(f, lA, lB, root_node -> children[i]);
    }
  }
}

int index_calculation(int i, int j, int L) {
  return (L - i) * (L - i + 1) / 2 + j;
}

void generate_part_0(FILE *f, char*variable) {
  fprintf(f, "         SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));\n");
  fprintf(f, "         SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));\n");
  fprintf(f, "         SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));\n");
  fprintf(f, "\n");
  fprintf(f, "         SCALAR_TYPE X_PC = SCALAR_SUB(x%s, xC);\n", variable);
  fprintf(f, "         SCALAR_TYPE Y_PC = SCALAR_SUB(y%s, yC);\n", variable);
  fprintf(f, "         SCALAR_TYPE Z_PC = SCALAR_SUB(z%s, zC);\n", variable);
  fprintf(f, "\n");
  fprintf(f, "         X_PC = SCALAR_MUL(X_PC, X_PC);\n");
  fprintf(f, "         X_PC = SCALAR_FMA(Y_PC, Y_PC, X_PC);\n");
  fprintf(f, "         X_PC = SCALAR_FMA(Z_PC, Z_PC, X_PC);\n");
  fprintf(f, "         SCALAR_TYPE TVAL = SCALAR_MUL(RHO, X_PC);\n\n");
}

void generate_part_1(FILE *f, int lA, int lB, struct node *root_node, char *variable) {  
  fprintf(f, "         SCALAR_TYPE ");
    
  for(int l = 1; l <= (lA + lB); ++l) {
    for(int k = 0; k < (lA + lB + 1) - l; ++k) {
      fprintf(f, "t%d%d, ", l, k);
    }
  }

  if((lA + lB) <= 1) {
    fprintf(f, "tx;\n\n");
  } else {
    fprintf(f, "tx, ty;\n\n");
  }

  for(int l = lA + lB - 1; l >= 0; --l) {
    fprintf(f, "         t0%d = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t0%d), TVAL_inv_e), SCALAR_SET1(%.20f));\n", l, l + 1, 2.0 / (1.0 * (2 * l + 1)));
  }
  fprintf(f, "\n");

  
  traverse_dfs_vrr(f, lA, lB, root_node);
}

void generate_diagonal_part_2(FILE *f, int lA, int type) {
  fprintf(f, "      double *Xik = (Xi + p_outer + p_inner);\n");
  fprintf(f, "      double *Gik = (Gi + p_outer + p_inner);\n");
  fprintf(f, "\n");

  if(type == 0) {
    fprintf(f, "      for(int c0 = 0; c0 <= %d; ++c0) {\n", lA);
    fprintf(f, "         for(int c1 = 0; c1 <= c0; ++c1) {\n");
    fprintf(f, "            int m = %d - c0;\n", lA);
    fprintf(f, "            int p = c1;\n");
    fprintf(f, "\n");
    fprintf(f, "            int idxB = (((%d - m) * (%d - m + 1)) >> 1) + p;\n", lA, lA);
    fprintf(f, "\n");
    fprintf(f, "            int mv, pv;\n");
    fprintf(f, "\n");

    fprintf(f, "            SCALAR_TYPE tx, wg, xik, gik;\n");
    
    int count = 0;
    for(int r0 = 0; r0 <= lA; ++r0) {
      for(int r1 = 0; r1 <= r0; ++r1) {
	int a = lA - r0;
	int c = r1;

	int idxA = index_calculation(a, c, lA);
	fprintf(f, "            mv = %d + m; pv = %d + p;\n", a, c);
	
	fprintf(f, "            tx  = SCALAR_LOAD((temp + (%d + (((%d - mv) * (%d - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));\n", (2 * lA * (2 * lA + 1) * (2 * lA + 2) - lA * (lA + 1) * (lA + 2)) / 6, 2 * lA, 2 * lA);
	fprintf(f, "            wg  = SCALAR_LOAD((weights + p_outer + p_inner));\n\n");
	fprintf(f, "            xik = SCALAR_LOAD((Xik + idxB * ldX));\n");
	fprintf(f, "            gik = SCALAR_LOAD((Gik + %d * ldG));\n\n", idxA);

	fprintf(f, "            tx = SCALAR_MUL(tx, wg);\n");
	fprintf(f, "            gik = SCALAR_FMA(tx, xik, gik);\n");
	fprintf(f, "            SCALAR_STORE((Gik + %d * ldG), gik);\n", idxA);

	count++;		
      }
    }
    fprintf(f, "         }\n");
    fprintf(f, "      }\n");
  } else if(type == 1) {
    fprintf(f, "      SCALAR_TYPE tx, wg, xik, gik;\n");
    
    for(int c0 = 0; c0 <= lA; ++c0) {
      for(int c1 = 0; c1 <= c0; ++c1) {
	int m = lA - c0;
	int p = c1;

	int idxB = index_calculation(m, p, lA);
	
	int count = 0;
	for(int r0 = 0; r0 <= lA; ++r0) {
	  for(int r1 = 0; r1 <= r0; ++r1) {
	    int a = lA - r0;
	    int c = r1;

	    int idxA = index_calculation(a, c, lA);

	    int idx = index_calculation(a + m - 0, c + p - 0, lA + lA - 0 - 0 - 0);

	    int offset = (2 * lA * (2 * lA + 1) * (2 * lA + 2) - lA * (lA + 1) * (lA + 2)) / 6;

	    fprintf(f, "      tx  = SCALAR_LOAD((temp + %d * blockDim.x + threadIdx.x));\n", offset + idx);
	    fprintf(f, "      wg  = SCALAR_LOAD((weights + p_outer + p_inner));\n\n");
	    fprintf(f, "      xik = SCALAR_LOAD((Xik + %d * ldX));\n", idxB);
	    fprintf(f, "      gik = SCALAR_LOAD((Gik + %d * ldG));\n\n", idxA);

	    fprintf(f, "      tx = SCALAR_MUL(tx, wg);\n");
	    fprintf(f, "      gik = SCALAR_FMA(tx, xik, gik);\n");
	    fprintf(f, "      SCALAR_STORE((Gik + %d * ldG), gik);\n", idxA);
      
	    count++;		
	  }
	}
      }
    }
  } else {
    printf("Type not defined\n");
  }
}

void generate_off_diagonal_part_2(FILE *f, int lA, int lB, int type) {
  fprintf(f, "      double *Xik = (Xi + p_outer + p_inner);\n");
  fprintf(f, "      double *Xjk = (Xj + p_outer + p_inner);\n");
  fprintf(f, "      double *Gik = (Gi + p_outer + p_inner);\n");
  fprintf(f, "      double *Gjk = (Gj + p_outer + p_inner);\n");
  fprintf(f, "\n");
  fprintf(f, "      SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));\n\n");
  
  if(type == 0) {
    fprintf(f, "      for(int c0 = 0; c0 <= %d; ++c0) {\n", lB);
    fprintf(f, "         for(int c1 = 0; c1 <= c0; ++c1) {\n");
    fprintf(f, "            int m = %d - c0;\n", lB);
    fprintf(f, "            int n = c0 - c1;\n");
    fprintf(f, "            int p = c1;\n");
    fprintf(f, "\n");
    fprintf(f, "            int idxB = (((%d - m) * (%d - m + 1)) >> 1) + p;\n", lB, lB);
    fprintf(f, "\n");
    fprintf(f, "            double X_ABp = 1.0, comb_m_i = 1.0;\n");
    fprintf(f, "            for(int i = 0; i <= m; ++i) {\n");
    fprintf(f, "               double rcp_i;\n");
    fprintf(f, "\n");
    fprintf(f, "               double Y_ABp = 1.0, comb_n_j = 1.0;\n");
    fprintf(f, "               for(int j = 0; j <= n; ++j) {\n");
    fprintf(f, "                  double rcp_j;\n");
    fprintf(f, "\n");
    fprintf(f, "                  double Z_ABp = 1.0, comb_p_k = 1.0;\n");
    fprintf(f, "                  for(int k = 0; k <= p; ++k) {\n");
    fprintf(f, "                     double rcp_k;\n");
    fprintf(f, "                     int mv, pv, Lv = %d - i - j - k;\n", lA + lB);
    fprintf(f, "\n");
    fprintf(f, "                     int offset = (Lv * (Lv + 1) * (Lv + 2) - %d) / 6;\n", lA * (lA + 1) * (lA + 2));
    fprintf(f, "                     double const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;\n");
    fprintf(f, "                     SCALAR_TYPE tx, ty, tz, tw;\n");
    fprintf(f, "                     SCALAR_TYPE const_value_w = SCALAR_MUL(const_value_v, const_value);\n\n");

    int count = 0;
    for(int r0 = 0; r0 <= lA; ++r0) {
      for(int r1 = 0; r1 <= r0; ++r1) {
	int a = lA - r0;
	int c = r1;

	int idxA = index_calculation(a, c, lA);
	fprintf(f, "                     mv = %d + m - i; pv = %d + p - k;\n", a, c);
	fprintf(f, "                     tx = SCALAR_LOAD((Xik + %d * ldX));\n", idxA);
	fprintf(f, "                     ty = SCALAR_LOAD((Xjk + idxB * ldX));\n");
	fprintf(f, "                     tz = SCALAR_LOAD((Gik + %d * ldG));\n", idxA);
	fprintf(f, "                     tw = SCALAR_LOAD((Gjk + idxB * ldG));\n");
	fprintf(f, "                     SCALAR_TYPE t%d = SCALAR_LOAD((temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * blockDim.x + threadIdx.x));\n", count);
	fprintf(f, "                     t%d = SCALAR_MUL(t%d, const_value_w);\n", count, count);
	fprintf(f, "                     tz = SCALAR_FMA(ty, t%d, tz);\n", count);
	fprintf(f, "                     tw = SCALAR_FMA(tx, t%d, tw);\n", count);
	fprintf(f, "                     SCALAR_STORE((Gik + %d * ldG), tz);\n", idxA);
	fprintf(f, "                     SCALAR_STORE((Gjk + idxB * ldG), tw);\n");
	count++;		
      }
    }
    fprintf(f, "\n");
    fprintf(f, "                     Z_ABp = SCALAR_MUL(Z_ABp, Z_AB);\n");
    fprintf(f, "                     rcp_k = SCALAR_RECIPROCAL(k + 1);\n");
    fprintf(f, "                     comb_p_k = SCALAR_MUL(comb_p_k, p - k);\n");
    fprintf(f, "                     comb_p_k = SCALAR_MUL(comb_p_k, rcp_k);\n");
    fprintf(f, "                  }\n");
    fprintf(f, "\n");
    fprintf(f, "                  Y_ABp = SCALAR_MUL(Y_ABp, Y_AB);\n");
    fprintf(f, "                  rcp_j = SCALAR_RECIPROCAL(j + 1);\n");
    fprintf(f, "                  comb_n_j = SCALAR_MUL(comb_n_j, n - j);\n");
    fprintf(f, "                  comb_n_j = SCALAR_MUL(comb_n_j, rcp_j);\n");
    fprintf(f, "               }\n");
    fprintf(f, "\n");
    fprintf(f, "               X_ABp = SCALAR_MUL(X_ABp, X_AB);\n");
    fprintf(f, "               rcp_i = SCALAR_RECIPROCAL(i + 1);\n");
    fprintf(f, "               comb_m_i = SCALAR_MUL(comb_m_i, m - i);\n");
    fprintf(f, "               comb_m_i = SCALAR_MUL(comb_m_i, rcp_i);\n");
    fprintf(f, "            }\n");
    fprintf(f, "         }\n");
    fprintf(f, "      }\n");
  } else if (type == 1) {
    fprintf(f, "      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;\n");
    fprintf(f, "      SCALAR_TYPE const_value_w;\n");

    int count = 0;
    fprintf(f, "      SCALAR_TYPE tx, ty, tz, tw, ");
    for(int r0 = 0; r0 <= lA - 1; ++r0) {
      for(int r1 = 0; r1 <= r0; ++r1) {
	fprintf(f, "t%d, ", count);
	count++;
      }
    }
    
    for(int r1 = 0; r1 <= lA - 1; ++r1) {
      fprintf(f, "t%d, ", count);
      count++;
    }

    fprintf(f, "t%d;\n", count);
    
    fprintf(f, "\n");

    for(int c0 = 0; c0 <= lB; ++c0) {
      for(int c1 = 0; c1 <= c0; ++c1) {
	int m = lB - c0;
	int n = c0 - c1;
	int p = c1;

	int idxB = index_calculation(m, p, lB);

	fprintf(f, "      X_ABp = 1.0; comb_m_i = 1.0;\n");
	for(int i = 0; i <= m; ++i) {
	  fprintf(f, "      Y_ABp = 1.0; comb_n_j = 1.0;\n");
	  for(int j = 0; j <= n; ++j) {
	    fprintf(f, "      Z_ABp = 1.0; comb_p_k = 1.0;\n");
	    for(int k = 0; k <= p; ++k) {
	      fprintf(f, "      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;\n");
	      fprintf(f, "      const_value_w = SCALAR_MUL(const_value_v, const_value);\n");

	      int count = 0;
	      for(int r0 = 0; r0 <= lA; ++r0) {
		for(int r1 = 0; r1 <= r0; ++r1) {
		  int a = lA - r0;
		  int c = r1;

		  int idxA = index_calculation(a, c, lA);

		  int idx = index_calculation(a + m - i, c + p - k, lA + lB - i - j - k);

		  int LAB = lA + lB - i - j - k;
		  int offset = (LAB * (LAB + 1) * (LAB + 2) - lA * (lA + 1) * (lA + 2)) / 6;
		  
		  fprintf(f, "      tx = SCALAR_LOAD((Xik + %d * ldX));\n", idxA);
		  fprintf(f, "      ty = SCALAR_LOAD((Xjk + %d * ldX));\n", idxB);
		  fprintf(f, "      tz = SCALAR_LOAD((Gik + %d * ldG));\n", idxA);
		  fprintf(f, "      tw = SCALAR_LOAD((Gjk + %d * ldG));\n", idxB);
		  fprintf(f, "      t%d = SCALAR_LOAD((temp + %d * blockDim.x + threadIdx.x));\n", count, offset + idx);
		  fprintf(f, "      t%d = SCALAR_MUL(t%d, const_value_w);\n", count, count);
		  fprintf(f, "      tz = SCALAR_FMA(ty, t%d, tz);\n", count);
		  fprintf(f, "      tw = SCALAR_FMA(tx, t%d, tw);\n", count);
		  fprintf(f, "      SCALAR_STORE((Gik + %d * ldG), tz);\n", idxA);
		  fprintf(f, "      SCALAR_STORE((Gjk + %d * ldG), tw);\n", idxB);
      
		  count++;		
		}
	      }
	      
	      if(k < p) {
		fprintf(f, "      Z_ABp = SCALAR_MUL(Z_ABp, Z_AB); comb_p_k = SCALAR_MUL(comb_p_k * %d, SCALAR_RECIPROCAL(%d));\n", p - k, k + 1);
	      }
	    }

	    if(j < n) {
	      fprintf(f, "      Y_ABp = SCALAR_MUL(Y_ABp, Y_AB); comb_n_j = SCALAR_MUL(comb_n_j * %d, SCALAR_RECIPROCAL(%d));\n", n - j, j + 1);
	    }
	  }

	  if(i < m) {
	    fprintf(f, "      X_ABp = SCALAR_MUL(X_ABp, X_AB); comb_m_i = SCALAR_MUL(comb_m_i * %d, SCALAR_RECIPROCAL(%d));\n", m - i, i + 1);
	  }
	}
      }
    }
  } else {
    printf("Type not defined\n");
  }  
}

void generate_diagonal_files(FILE *f, int lA, int size, struct node *root_node, int type) {
  fprintf(f, "#include <math.h>\n");
  fprintf(f, "#include \"../include/chebyshev_boys_computation.hpp\"\n");
  fprintf(f, "#include \"../include/integral_data_types.hpp\"\n");
  fprintf(f, "#include \"config_obara_saika.hpp\"\n");
  fprintf(f, "#include \"integral_%d.hu\"\n", lA);
  fprintf(f, "\n");
  fprintf(f, "#define PI 3.14159265358979323846\n");
  fprintf(f, "\n");
  fprintf(f, "#define MIN(a,b)			\\\n"); 
  fprintf(f, "  ({ __typeof__ (a) _a = (a);	        \\\n");
  fprintf(f, "  __typeof__ (b) _b = (b);		\\\n");
  fprintf(f, "  _a < _b ? _a : _b; })\n");
  fprintf(f, "\n");
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "__global__ void integral_%d(size_t npts,\n", lA);
  fprintf(f, "                          point rA,\n");
  fprintf(f, "                          point rB,\n");
  fprintf(f, "                          int nprim_pairs,\n");
  fprintf(f, "                          prim_pair *prim_pairs,\n");  
  fprintf(f, "                          double *_points,\n");
  fprintf(f, "                          double *Xi,\n");
  fprintf(f, "                          int ldX,\n");
  fprintf(f, "                          double *Gi,\n");
  fprintf(f, "                          int ldG, \n");
  fprintf(f, "                          double *weights,\n");
  fprintf(f, "                          double *boys_table) {\n");	 

  int partial_size = 0;
  for(int i = 0; i < lA; ++i) {
    partial_size += (i + 1) * (i + 2) / 2;
  }
  
  //fprintf(f, "   __shared__ double temp[%d * blockDim.x];\n",  size - partial_size);
  fprintf(f, "   __shared__ double *temp;\n");
  
  char variable[1024];
  sprintf(variable, "A");
  
  fprintf(f, "   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {\n");
  fprintf(f, "      double *_point_outer = (_points + p_outer);\n\n");
  fprintf(f, "      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);\n\n");
  fprintf(f, "      double xA = rA.x;\n");
  fprintf(f, "      double yA = rA.y;\n");
  fprintf(f, "      double zA = rA.z;\n");
  fprintf(f, "\n");
  fprintf(f, "      for(int i = 0; i < %d; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());\n", size - partial_size);
  fprintf(f, "\n");
  fprintf(f, "      for(int ij = 0; ij < nprim_pairs; ++ij) {\n");
  fprintf(f, "         double RHO = prim_pairs[ij].gamma;\n");
  if(lA > 0) {
    fprintf(f, "         double RHO_INV = prim_pairs[ij].gamma_inv;\n");
  }
  fprintf(f, "\n");
  if(lA != 0) {
    fprintf(f, "         constexpr double X_PA = 0.0;\n");
    fprintf(f, "         constexpr double Y_PA = 0.0;\n");
    fprintf(f, "         constexpr double Z_PA = 0.0;\n");
    fprintf(f, "\n");
  }
  fprintf(f, "         double eval = prim_pairs[ij].K_coeff_prod;\n");
  fprintf(f, "\n");
  fprintf(f, "         // Evaluate T Values\n");
  
  generate_part_0(f, variable);

  fprintf(f, "         SCALAR_TYPE ");
  for(int l = 0; l < 2 * lA; ++l) {
    fprintf(f, "t0%d, ", l);
  }
  fprintf(f, "t0%d, TVAL_inv_e;\n\n", 2 * lA);

  fprintf(f, "         // Evaluate Boys function\n");
  fprintf(f, "         boys_element<%d>(&TVAL, &TVAL_inv_e, &t0%d, boys_table);\n", 2 * lA, 2 * lA);
  fprintf(f, "\n");
  fprintf(f, "         // Evaluate VRR Buffer\n");

  generate_part_1(f, lA, lA, root_node, variable);
  
  fprintf(f, "      }\n");
  fprintf(f, "\n");  

  generate_diagonal_part_2(f, lA, type);
  
  fprintf(f, "   }\n");
  fprintf(f, "}\n");
  fprintf(f, "}\n");
}

void generate_off_diagonal_files(FILE *f, int lA, int lB, int size, struct node *root_node, int type) {
  fprintf(f, "#include <math.h>\n");
  fprintf(f, "#include \"../include/chebyshev_boys_computation.hpp\"\n");
  fprintf(f, "#include \"../include/integral_data_types.hpp\"\n");
  fprintf(f, "#include \"config_obara_saika.hpp\"\n");
  fprintf(f, "#include \"integral_%d_%d.hu\"\n", lA, lB);
  fprintf(f, "\n");
  fprintf(f, "#define PI 3.14159265358979323846\n");
  fprintf(f, "\n");
  fprintf(f, "#define MIN(a,b)			\\\n"); 
  fprintf(f, "  ({ __typeof__ (a) _a = (a);	        \\\n");
  fprintf(f, "  __typeof__ (b) _b = (b);		\\\n");
  fprintf(f, "  _a < _b ? _a : _b; })\n");
  fprintf(f, "\n");
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "__global__ void integral_%d_%d(size_t npts,\n", lA, lB);
  fprintf(f, "                             point rA,\n");
  fprintf(f, "                             point rB,\n");
  fprintf(f, "                             int nprim_pairs,\n");
  fprintf(f, "                             prim_pair *prim_pairs,\n");
  fprintf(f, "                             double *_points,\n");
  fprintf(f, "                             double *Xi,\n");
  fprintf(f, "                             double *Xj,\n");
  fprintf(f, "                             int ldX,\n");
  fprintf(f, "                             double *Gi,\n");
  fprintf(f, "                             double *Gj,\n");
  fprintf(f, "                             int ldG, \n");
  fprintf(f, "                             double *weights, \n");
  fprintf(f, "                             double *boys_table) {\n");	 

  int partial_size = 0;
  for(int i = 0; i < lA; ++i) {
    partial_size += (i + 1) * (i + 2) / 2;
  }
  
  //fprintf(f, "   __shared__ double temp[%d * blockDim.x];\n",  size - partial_size);
  fprintf(f, "   __shared__ double *temp;\n");

  char variable[1024];
  sprintf(variable, "P");

  fprintf(f, "   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {\n");
  fprintf(f, "      double *_point_outer = (_points + p_outer);\n\n");
  fprintf(f, "      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);\n\n");
  if(lB != 0) {
    fprintf(f, "      double X_AB = rA.x - rB.x;\n");
    fprintf(f, "      double Y_AB = rA.y - rB.y;\n");
    fprintf(f, "      double Z_AB = rA.z - rB.z;\n");
    fprintf(f, "\n");
  }
  fprintf(f, "      for(int i = 0; i < %d; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());\n", size - partial_size);
  fprintf(f, "\n");
  fprintf(f, "      for(int ij = 0; ij < nprim_pairs; ++ij) {\n");
  fprintf(f, "         double RHO = prim_pairs[ij].gamma;\n");
  if(lA + lB > 1) {
    fprintf(f, "         double RHO_INV = prim_pairs[ij].gamma_inv;\n");
  }
  if(lA != 0) {
    fprintf(f, "         double X_PA = prim_pairs[ij].PA.x;\n");
    fprintf(f, "         double Y_PA = prim_pairs[ij].PA.y;\n");
    fprintf(f, "         double Z_PA = prim_pairs[ij].PA.z;\n");
  }
  fprintf(f, "\n");
  fprintf(f, "         double xP = prim_pairs[ij].P.x;\n");
  fprintf(f, "         double yP = prim_pairs[ij].P.y;\n");
  fprintf(f, "         double zP = prim_pairs[ij].P.z;\n");
  fprintf(f, "\n");
  fprintf(f, "         double eval = prim_pairs[ij].K_coeff_prod;\n");
  fprintf(f, "\n");
  fprintf(f, "         // Evaluate T Values\n");

  generate_part_0(f, variable);

  fprintf(f, "         SCALAR_TYPE ");
  for(int l = 0; l < lA+lB; ++l) {
    fprintf(f, "t0%d, ", l);
  }
  fprintf(f, "t0%d, TVAL_inv_e;\n\n", (lA + lB));
  
  fprintf(f, "         // Evaluate Boys function\n");
  fprintf(f, "         boys_element<%d>(&TVAL, &TVAL_inv_e, &t0%d, boys_table);\n", lA + lB, lA + lB);
  fprintf(f, "\n");
  fprintf(f, "         // Evaluate VRR Buffer\n");

  generate_part_1(f, lA, lB, root_node, variable);
  
  fprintf(f, "      }\n");
  fprintf(f, "\n");

  generate_off_diagonal_part_2(f, lA, lB, type);

  fprintf(f, "   }\n");
  fprintf(f, "}\n");
  fprintf(f, "}\n");
}

void generate_diagonal_header_files(int lA) {
  char filename[512];
      
  sprintf(filename, "integral_%d.hu", lA);
      
  FILE *f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_%d\n", lA);
  fprintf(f, "#define __MY_INTEGRAL_%d\n", lA);
  fprintf(f, "\n");
  fprintf(f, "#include \"../include/integral_data_types.hpp\"\n");
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "__global__ void integral_%d(size_t npts,\n", lA);
  fprintf(f, "                          point rA,\n");
  fprintf(f, "                          point rB,\n");
  fprintf(f, "                          int nprim_pairs,\n");
  fprintf(f, "                          prim_pair *prim_pairs,\n");
  fprintf(f, "                          double *points,\n");
  fprintf(f, "                          double *Xi,\n");
  fprintf(f, "                          int ldX,\n");	 
  fprintf(f, "                          double *Gi,\n");
  fprintf(f, "                          int ldG, \n");
  fprintf(f, "                          double *weights, \n");
  fprintf(f, "                          double *boys_table);\n");	 
  fprintf(f, "}\n");
  fprintf(f, "#endif\n");
  
  fclose(f);
}

void generate_off_diagonal_header_files(int lA, int lB) {
  char filename[512];
      
  sprintf(filename, "integral_%d_%d.hu", lA, lB);
      
  FILE *f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_%d_%d\n", lA, lB);
  fprintf(f, "#define __MY_INTEGRAL_%d_%d\n", lA, lB);
  fprintf(f, "\n");
  fprintf(f, "#include \"../include/integral_data_types.hpp\"\n");
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "__global__ void integral_%d_%d(size_t npts,\n", lA, lB);
  fprintf(f, "                             point rA,\n");
  fprintf(f, "                             point rB,\n");
  fprintf(f, "                             int nprim_pairs,\n");
  fprintf(f, "                             prim_pair *prim_pairs,\n");
  fprintf(f, "                             double *points,\n");
  fprintf(f, "                             double *Xi,\n");
  fprintf(f, "                             double *Xj,\n");
  fprintf(f, "                             int ldX,\n");	 
  fprintf(f, "                             double *Gi,\n");
  fprintf(f, "                             double *Gj,\n");
  fprintf(f, "                             int ldG, \n");
  fprintf(f, "                             double *weights,\n");
  fprintf(f, "                             double *boys_table);\n");	 
  fprintf(f, "}\n");
  fprintf(f, "#endif\n");
  
  fclose(f);
}

void generate_main_files(int lA) {
  char filename[512];

  FILE *f;
  
  sprintf(filename, "obara_saika_integrals.hpp");
      
  f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_OBARA_SAIKA\n");
  fprintf(f, "#define __MY_INTEGRAL_OBARA_SAIKA\n");
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "void generate_shell_pair( const shells& A, const shells& B, prim_pair *prim_pairs);\n");
  fprintf(f, "void compute_integral_shell_pair(size_t npts,\n");
  fprintf(f, "                             int is_diag,\n");
  fprintf(f, "                             int lA,\n");
  fprintf(f, "                             int lB,\n");
  fprintf(f, "                             point rA,\n");
  fprintf(f, "                             point rB,\n");
  fprintf(f, "                             int nprim_pairs,\n");
  fprintf(f, "                             prim_pair *prim_pairs,\n");
  fprintf(f, "                             double *points,\n");
  fprintf(f, "                             double *Xi,\n");
  fprintf(f, "                             double *Xj,\n");
  fprintf(f, "                             int ldX,\n");	 
  fprintf(f, "                             double *Gi,\n");
  fprintf(f, "                             double *Gj,\n");
  fprintf(f, "                             int ldG, \n");
  fprintf(f, "                             double *weights,\n");
  fprintf(f, "                             double *boys_table);\n");
  fprintf(f, " }\n");
  fprintf(f, "#endif\n");
  
  fclose(f);  

  sprintf(filename, "obara_saika_integrals.cu");
      
  f = fopen(filename, "w");

  fprintf(f, "#include <stdio.h>\n");
  fprintf(f, "#include <stdlib.h>\n");
  fprintf(f, "#include \"../include/integral_data_types.hpp\"\n");
  fprintf(f, "#include \"../include/obara_saika_integrals.hpp\"\n");
  for(int i = 0; i <= lA; ++i) {
    fprintf(f, "#include \"integral_%d.hu\"\n", i);
  }

  for(int i = 0; i <= lA; ++i) {
    for(int j = 0; j <= i; ++j) {
      fprintf(f, "#include \"integral_%d_%d.hu\"\n", i, j);
    }
  }
  fprintf(f, "namespace XGPU {\n");
  fprintf(f, "\nvoid generate_shell_pair( const shells& A, const shells& B, prim_pair *prim_pairs) {\n");
  fprintf(f, "   // L Values\n");
  fprintf(f, "   int lA = A.L;\n");
  fprintf(f, "   int lB = B.L;\n\n");

  fprintf(f, "   const auto xA = A.origin.x;\n");
  fprintf(f, "   const auto yA = A.origin.y;\n");
  fprintf(f, "   const auto zA = A.origin.z;\n\n");

  fprintf(f, "   const auto xB = B.origin.x;\n");
  fprintf(f, "   const auto yB = B.origin.y;\n");
  fprintf(f, "   const auto zB = B.origin.z;\n\n");

  fprintf(f, "   double rABx = xA - xB;\n");
  fprintf(f, "   double rABy = yA - yB;\n");
  fprintf(f, "   double rABz = zA - zB;\n\n");

  fprintf(f, "   const double dAB = rABx*rABx + rABy*rABy + rABz*rABz;\n\n");

  fprintf(f, "   const int nprim_A = A.m;\n");
  fprintf(f, "   const int nprim_B = B.m;\n");
  fprintf(f, "   const int np = nprim_A * nprim_B;\n\n");

  fprintf(f, "   for(int i = 0, ij = 0; i < nprim_A; ++i       )\n");
  fprintf(f, "   for(int j = 0        ; j < nprim_B; ++j, ++ij ) {\n");
  fprintf(f, "      auto& pair = prim_pairs[ij];\n");
  fprintf(f, "      const auto alpha_A = A.coeff[i].alpha;\n");
  fprintf(f, "      const auto alpha_B = B.coeff[j].alpha;\n\n");

  fprintf(f, "      pair.gamma = alpha_A + alpha_B;\n");
  fprintf(f, "      pair.gamma_inv = 1. / pair.gamma;\n\n");

  fprintf(f, "      pair.P.x = (alpha_A * xA + alpha_B * xB) * pair.gamma_inv;\n");
  fprintf(f, "      pair.P.y = (alpha_A * yA + alpha_B * yB) * pair.gamma_inv;\n");
  fprintf(f, "      pair.P.z = (alpha_A * zA + alpha_B * zB) * pair.gamma_inv;\n\n");

  fprintf(f, "      pair.PA.x = pair.P.x - xA;\n");
  fprintf(f, "      pair.PA.y = pair.P.y - yA;\n");
  fprintf(f, "      pair.PA.z = pair.P.z - zA;\n\n");

  fprintf(f, "      pair.PB.x = pair.P.x - xB;\n");
  fprintf(f, "      pair.PB.y = pair.P.y - yB;\n");
  fprintf(f, "      pair.PB.z = pair.P.z - zB;\n\n");

  fprintf(f, "      pair.K_coeff_prod = 2 * M_PI * A.coeff[i].coeff * B.coeff[j].coeff * pair.gamma_inv * std::exp( - alpha_A * alpha_B * dAB * pair.gamma_inv );\n");
  fprintf(f, "   }\n\n");
  fprintf(f, "}\n");
  
  fprintf(f, "\n");
  fprintf(f, "void compute_integral_shell_pair(size_t npts,\n");
  fprintf(f, "                  int is_diag,\n");
  fprintf(f, "                  int lA,\n");
  fprintf(f, "                  int lB,\n");
  fprintf(f, "                  point rA,\n");
  fprintf(f, "                  point rB,\n");
  fprintf(f, "                  int nprim_pairs,\n");
  fprintf(f, "                  prim_pair *prim_pairs,\n");
  fprintf(f, "                  double *points,\n");
  fprintf(f, "                  double *Xi,\n");
  fprintf(f, "                  double *Xj,\n");
  fprintf(f, "                  int ldX,\n");	 
  fprintf(f, "                  double *Gi,\n");
  fprintf(f, "                  double *Gj,\n");
  fprintf(f, "                  int ldG, \n");
  fprintf(f, "                  double *weights,\n");
  fprintf(f, "                  double *boys_table) {\n");	 

  int size = 0;
  int partial_size = 0;
  
  fprintf(f, "   if (is_diag) {\n");
  size = 0;
  for(int l = 0; l < (0 + 0 + 1); ++l) {
    size += (l + 1) * (l + 2) / 2;
  }

  partial_size = 0;
  for(int l = 0; l < 0; ++l) {
    partial_size += (l + 1) * (l + 2) / 2;
  }
  
  fprintf(f, "      if(lA == %d) {\n", 0);
  fprintf(f, "         integral_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", 0, size - partial_size);
  fprintf(f, "                                rA,\n");
  fprintf(f, "                                rB,\n");
  fprintf(f, "                                nprim_pairs,\n");
  fprintf(f, "                                prim_pairs,\n");
  fprintf(f, "                                points,\n");
  fprintf(f, "                                Xi,\n");
  fprintf(f, "                                ldX,\n");
  fprintf(f, "                                Gi,\n");
  fprintf(f, "                                ldG, \n");
  fprintf(f, "                                weights, \n");
  fprintf(f, "                                boys_table);\n");
  
  fprintf(f, "      } else ");

  for(int i = 1; i <= lA; ++i) {
    size = 0;
    for(int l = 0; l < (i + i + 1); ++l) {
      size += (l + 1) * (l + 2) / 2;
    }

    partial_size = 0;
    for(int l = 0; l < i; ++l) {
      partial_size += (l + 1) * (l + 2) / 2;
    }

    fprintf(f, "if(lA == %d) {\n", i);
    fprintf(f, "        integral_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", 0, size - partial_size);
    fprintf(f, "                               rA,\n");
    fprintf(f, "                               rB,\n");
    fprintf(f, "                               nprim_pairs,\n");
    fprintf(f, "                               prim_pairs,\n");
    fprintf(f, "                               points,\n");
    fprintf(f, "                               Xi,\n");
    fprintf(f, "                               ldX,\n");
    fprintf(f, "                               Gi,\n");
    fprintf(f, "                               ldG, \n");
    fprintf(f, "                               weights, \n");
    fprintf(f, "                               boys_table);\n");	 
    fprintf(f, "      } else ");
  }

  fprintf(f, "{\n");
  fprintf(f, "         printf(\"Type not defined!\\n\");\n");
  fprintf(f, "      }\n");  
  fprintf(f, "   } else {\n");

  size = 0;
  for(int l = 0; l < (0 + 0 + 1); ++l) {
    size += (l + 1) * (l + 2) / 2;
  }

  partial_size = 0;
  for(int l = 0; l < 0; ++l) {
    partial_size += (l + 1) * (l + 2) / 2;
  }
  
  fprintf(f, "      if((lA == %d) && (lB == %d)) {\n", 0, 0);
  fprintf(f, "         integral_%d_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", 0, 0, size - partial_size);
  fprintf(f, "                                  rA,\n");
  fprintf(f, "                                  rB,\n");
  fprintf(f, "                                  nprim_pairs,\n");
  fprintf(f, "                                  prim_pairs,\n");
  fprintf(f, "                                  points,\n");
  fprintf(f, "                                  Xi,\n");
  fprintf(f, "                                  Xj,\n");
  fprintf(f, "                                  ldX,\n");
  fprintf(f, "                                  Gi,\n");
  fprintf(f, "                                  Gj,\n");
  fprintf(f, "                                  ldG, \n");
  fprintf(f, "                                  weights,\n");
  fprintf(f, "                                  boys_table);\n");	 
  fprintf(f, "      } else ");

  for(int i = 1; i <= lA; ++i) {
    for(int j = 0; j < i; ++j) {
      size = 0;
      for(int l = 0; l < (i + j + 1); ++l) {
	size += (l + 1) * (l + 2) / 2;
      }

      partial_size = 0;
      for(int l = 0; l < i; ++l) {
	partial_size += (l + 1) * (l + 2) / 2;
      }
      
      fprintf(f, "if((lA == %d) && (lB == %d)) {\n", i, j);
      fprintf(f, "         integral_%d_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", i, j, size - partial_size);
      fprintf(f, "                                  rA,\n");
      fprintf(f, "                                  rB,\n");
      fprintf(f, "                                  nprim_pairs,\n");
      fprintf(f, "                                  prim_pairs,\n");
      fprintf(f, "                                  points,\n");
      fprintf(f, "                                  Xi,\n");
      fprintf(f, "                                  Xj,\n");
      fprintf(f, "                                  ldX,\n");
      fprintf(f, "                                  Gi,\n");
      fprintf(f, "                                  Gj,\n");
      fprintf(f, "                                  ldG, \n");
      fprintf(f, "                                  weights,\n");
      fprintf(f, "                                  boys_table);\n");	 
      fprintf(f, "      } else if((lA == %d) && (lB == %d)) {\n", j, i);
      fprintf(f, "         integral_%d_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", i, j, size - partial_size);
      fprintf(f, "                                  rB,\n");
      fprintf(f, "                                  rA,\n");
      fprintf(f, "                                  nprim_pairs,\n");
      fprintf(f, "                                  prim_pairs,\n");
      fprintf(f, "                                  points,\n");
      fprintf(f, "                                  Xj,\n");
      fprintf(f, "                                  Xi,\n");
      fprintf(f, "                                  ldX,\n");
      fprintf(f, "                                  Gj,\n");
      fprintf(f, "                                  Gi,\n");
      fprintf(f, "                                  ldG, \n");
      fprintf(f, "                                  weights, \n");
      fprintf(f, "                                  boys_table);\n");	 
      fprintf(f, "      } else ");
    }

    size = 0;
    for(int l = 0; l < (i + i + 1); ++l) {
      size += (l + 1) * (l + 2) / 2;
    }

    partial_size = 0;
    for(int l = 0; l < i; ++l) {
      partial_size += (l + 1) * (l + 2) / 2;
    }
    
    fprintf(f, "if((lA == %d) && (lB == %d)) {\n", i, i);
    fprintf(f, "        integral_%d_%d<<<320, 128, 128 * %d * sizeof(double)>>>(npts,\n", i, i, size - partial_size);
    fprintf(f, "                                 rA,\n");
    fprintf(f, "                                 rB,\n");
    fprintf(f, "                                 nprim_pairs,\n");
    fprintf(f, "                                 prim_pairs,\n");
    fprintf(f, "                                 points,\n");
    fprintf(f, "                                 Xi,\n");
    fprintf(f, "                                 Xj,\n");
    fprintf(f, "                                 ldX,\n");
    fprintf(f, "                                 Gi,\n");
    fprintf(f, "                                 Gj,\n");
    fprintf(f, "                                 ldG, \n");
    fprintf(f, "                                 weights,\n");
    fprintf(f, "                                 boys_table);\n");	 
    fprintf(f, "      } else ");
  }

  fprintf(f, "{\n");
  fprintf(f, "         printf(\"Type not defined!\\n\");\n");
  fprintf(f, "      }\n");
  fprintf(f, "   }\n\n");  
  fprintf(f, "}\n");
  fprintf(f, "}\n");
  
  fclose(f);  
}

int main(int argc, char **argv) {
  int lA = atoi(argv[1]);
  int tV = atoi(argv[2]);

  generate_main_files(lA);
  
  for(int i = 0; i <= lA; ++i) {
    for(int j = 0; j <= i; ++j) {
      int size = 0;
      for(int l = 0; l < (i + j + 1); ++l) {
	size += (l + 1) * (l + 2) / 2;
      }
  
      struct node *node_list = (struct node *) malloc(size * sizeof(struct node));

      for(int i = 0; i < size; ++i) {
	node_list[i].iA = 0;
	node_list[i].jA = 0;
	node_list[i].kA = 0;
    
	node_list[i].iB = 0;
	node_list[i].jB = 0;
	node_list[i].kB = 0;
    
	node_list[i].level = 0;
	node_list[i].vars = 0;

	node_list[i].valid = 0;
	node_list[i].offset = 0;
    
	node_list[i].nr_children = 0;;
      }

      if(i == j) {
	generate_diagonal_header_files(i);
      }
      
      generate_off_diagonal_header_files(i, j);
      
      int type = ((i + j) <= tV) ? 1 : 0;
      
      // initialization part
      initialize_tree_structure(type, i, j, size, node_list);
      
      // vrr construction
      if(i == j) {
	char filename[512];
      
	sprintf(filename, "integral_%d.cu", i);
      
	FILE *f = fopen(filename, "w");

	generate_diagonal_files(f, i, size, node_list, type);

	fclose(f);
      }

      char filename[512];
      
      sprintf(filename, "integral_%d_%d.cu", i, j);
      
      FILE *f = fopen(filename, "w");

      generate_off_diagonal_files(f, i, j, size, node_list, type);

      fclose(f);
      
      free(node_list);
    }
  }
  
  return 0;
}
