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
	fprintf(f, "            t%d%d = eval * boys_function(%d, tval);\n", root_node -> level, v, v);
      }
    } else if (root_node -> level == 1) {
      for(int v = 0; v < root_node -> vars; ++v) {
	fprintf(f, "            t%d%d = %s * t%d%d - %s * t%d%d;\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v, root_node -> var_pc, root_node -> level - 1, v + 1);
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
	  fprintf(f, "            t%d%d = %s * t%d%d - %s * t%d%d;\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v, root_node -> var_pc, root_node -> level - 1, v + 1);
	}
      } else {
	for(int v = 0; v < root_node -> vars; ++v) {
	  fprintf(f, "            t%d%d = %s * t%d%d - %s * t%d%d + 0.5 * RHO_INV * %d * (t%d%d - t%d%d);\n", root_node -> level, v, root_node -> var_pa, root_node -> level - 1, v, root_node -> var_pc, root_node -> level - 1, v + 1, iteration, root_node -> level - 2, v, root_node ->level - 2, v + 1);
	}
      }
    }

    if(root_node -> valid) {
      fprintf(f, "            *(temp + %d * NPTS_LOCAL + p_inner) += t%d%d;\n", root_node -> offset, root_node -> level, 0);
    }
    
    for(int i = 0; i < root_node -> nr_children; ++i) {
      traverse_dfs_vrr(f, lA, lB, root_node -> children[i]);
    }
  }
}

int index_calculation(int i, int j, int L) {
  return (L - i) * (L - i + 1) / 2 + j;
}

void generate_diagonal_files(FILE *f, int lA, int size, struct node *root_node, int type) {
  fprintf(f, "#include <math.h>\n");
  fprintf(f, "#include \"boys_computation.h\"\n");
  fprintf(f, "#include \"integral_data_types.h\"\n");
  fprintf(f, "\n");
  fprintf(f, "#define PI 3.14159265358979323846\n");
  fprintf(f, "\n");
  fprintf(f, "#define MIN(a,b)			\\\n"); 
  fprintf(f, "  ({ __typeof__ (a) _a = (a);	        \\\n");
  fprintf(f, "  __typeof__ (b) _b = (b);		\\\n");
  fprintf(f, "  _a < _b ? _a : _b; })\n");
  fprintf(f, "\n");
  fprintf(f, "void integral_%d(size_t npts,\n", lA);
//fprintf(f, "               shells shellA,\n");
  fprintf(f, "               shell_pair shpair,\n");
  fprintf(f, "               point *_points,\n");
  fprintf(f, "               double *Xi,\n");
  fprintf(f, "               int stX,\n");
  fprintf(f, "               int ldX,\n");
  fprintf(f, "               double *Gi,\n");
  fprintf(f, "               int stG, \n");
  fprintf(f, "               int ldG, \n");
  fprintf(f, "               double *weights) {\n");	 

  int partial_size = 0;
  for(int i = 0; i < lA; ++i) {
    partial_size += (i + 1) * (i + 2) / 2;
  }

  fprintf(f, "   double temp[%d * NPTS_LOCAL];\n\n", size - partial_size);
  fprintf(f, "   for(int i = 0; i < %d * NPTS_LOCAL; ++i) {\n", size - partial_size);
  fprintf(f, "      temp[i] = 0.0;\n");
  fprintf(f, "   }\n\n");
  
  fprintf(f, "   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {\n");
  fprintf(f, "      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);\n");
  fprintf(f, "      point *_point_outer = (_points + p_outer);\n\n");

//fprintf(f, "      double xA = shellA.origin.x;\n");
//fprintf(f, "      double yA = shellA.origin.y;\n");
//fprintf(f, "      double zA = shellA.origin.z;\n");
  fprintf(f, "      double xA = shpair.rA.x;\n");
  fprintf(f, "      double yA = shpair.rA.y;\n");
  fprintf(f, "      double zA = shpair.rA.z;\n");
  fprintf(f, "\n");
//fprintf(f, "      double beta_in = 0.0;\n");
//fprintf(f, "      for(int i = 0; i < shellA.m; ++i) {\n");
//fprintf(f, "         for(int j = 0; j < shellA.m; ++j) {\n");
  fprintf(f, "      for( int ij = 0; ij < shpair.nprim_pair; ++ij ) {\n");
//fprintf(f, "         double aA = shellA.coeff[i].alpha;\n");
//fprintf(f, "         double cA = shellA.coeff[i].coeff;\n");
//fprintf(f, "\n");
//fprintf(f, "         double aB = shellA.coeff[j].alpha;\n");
//fprintf(f, "         double cB = shellA.coeff[j].coeff;\n");
//fprintf(f, "\n");
//fprintf(f, "         double RHO = aA + aB;\n");
  fprintf(f, "         double RHO = shpair.prim_pairs[ij].gamma;\n");
  fprintf(f, "         double RHO_INV = 1.0 / RHO;\n");
  fprintf(f, "\n");  
  fprintf(f, "         constexpr double X_PA = 0.0;\n");
  fprintf(f, "         constexpr double Y_PA = 0.0;\n");
  fprintf(f, "         constexpr double Z_PA = 0.0;\n");
  fprintf(f, "\n");
  fprintf(f, "         double eval = shpair.prim_pairs[ij].coeff_prod * 2 * PI * RHO_INV;\n");
  fprintf(f, "\n");
  fprintf(f, "         for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {\n");
  fprintf(f, "            point C = *(_point_outer + p_inner);\n");
  fprintf(f, "\n");  
  fprintf(f, "            double xC = C.x;\n");
  fprintf(f, "            double yC = C.y;\n");
  fprintf(f, "            double zC = C.z;\n");
  fprintf(f, "\n");
  fprintf(f, "            double X_PC = (xA - xC);\n");
  fprintf(f, "            double Y_PC = (yA - yC);\n");
  fprintf(f, "            double Z_PC = (zA - zC);\n");
  fprintf(f, "\n");
  fprintf(f, "            double ");
  for(int l = 0; l < (lA + lA); ++l) {
    for(int k = 0; k < (lA + lA + 1) - l; ++k) {
      fprintf(f, "t%d%d, ", l, k);
    }
  }
  fprintf(f, "t%d%d;\n", (lA + lA), 0);
  fprintf(f, "\n");
  fprintf(f, "            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);\n");
  fprintf(f, "\n");
  traverse_dfs_vrr(f, lA, lA, root_node);
//fprintf(f, "            beta_in = 1.0;\n");
  fprintf(f, "         }\n");
  fprintf(f, "      }\n");
  fprintf(f, "\n");
  fprintf(f, "      for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {;\n");
  fprintf(f, "         double *Xik = (Xi + (NPTS_LOCAL * p_outer + p_inner) * stX);\n");
  fprintf(f, "         double *Gik = (Gi + (NPTS_LOCAL * p_outer + p_inner) * stG);\n");
  fprintf(f, "\n");

  if(type == 0) {
    fprintf(f, "         for(int c0 = 0; c0 <= %d; ++c0) {\n", lA);
    fprintf(f, "            for(int c1 = 0; c1 <= c0; ++c1) {\n");
    fprintf(f, "               int m = %d - c0;\n", lA);
    fprintf(f, "               int p = c1;\n");
    fprintf(f, "\n");
    fprintf(f, "               int idxB = (((%d - m) * (%d - m + 1)) >> 1) + p;\n", lA, lA);
    fprintf(f, "\n");
    fprintf(f, "               int mv, pv;\n");
    fprintf(f, "\n");
    
    int count = 0;
    for(int r0 = 0; r0 <= lA; ++r0) {
      for(int r1 = 0; r1 <= r0; ++r1) {
	int a = lA - r0;
	int c = r1;

	int idxA = index_calculation(a, c, lA);
	fprintf(f, "               mv = %d + m; pv = %d + p;\n", a, c);
	fprintf(f, "               *(Gik + %d * ldG) += *(Xik + idxB * ldX) * (*(temp + (%d + (((%d - mv) * (%d - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner)) * (*(weights + (NPTS_LOCAL * p_outer + p_inner)));\n", idxA, (2 * lA * (2 * lA + 1) * (2 * lA + 2) - lA * (lA + 1) * (lA + 2)) / 6, 2 * lA, 2 * lA);
      
	//if (idxA != ((lA + 1) * (lA + 2) / 2 - 1)) fprintf(f, "\n");
	count++;		
      }
    }
    fprintf(f, "            }\n");
    fprintf(f, "         }\n");
    fprintf(f, "      }\n");
  } else if(type == 1) {   
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
		  
	    fprintf(f, "         *(Gik + %d * ldG) += *(Xik + %d * ldX) * (*(temp + %d * NPTS_LOCAL + p_inner)) * (*(weights + (NPTS_LOCAL * p_outer + p_inner)));\n", idxA, idxB, offset + idx);
      
	    count++;		
	  }
	}

	//if(idxB != ((lA + 1) * (lA + 2) / 2 - 1)) fprintf(f, "\n");
      }
    }
    fprintf(f, "      }\n");
  } else {
    fprintf(f, "Type not defined\n");
  }  

  fprintf(f, "   }\n");
  fprintf(f, "}\n");
}

void generate_off_diagonal_files(FILE *f, int lA, int lB, int size, struct node *root_node, int type) {
  fprintf(f, "#include <math.h>\n");
  fprintf(f, "#include \"boys_computation.h\"\n");
  fprintf(f, "#include \"integral_data_types.h\"\n");
  fprintf(f, "\n");
  fprintf(f, "#define PI 3.14159265358979323846\n");
  fprintf(f, "\n");
  fprintf(f, "#define MIN(a,b)			\\\n"); 
  fprintf(f, "  ({ __typeof__ (a) _a = (a);	        \\\n");
  fprintf(f, "  __typeof__ (b) _b = (b);		\\\n");
  fprintf(f, "  _a < _b ? _a : _b; })\n");
  fprintf(f, "\n");
  fprintf(f, "void integral_%d_%d(size_t npts,\n", lA, lB);
//fprintf(f, "                  shells shellA,\n");
//fprintf(f, "                  shells shellB,\n");
  fprintf(f, "                  shell_pair shpair,\n");
  fprintf(f, "                  point *_points,\n");
  fprintf(f, "                  double *Xi,\n");
  fprintf(f, "                  double *Xj,\n");
  fprintf(f, "                  int stX,\n");
  fprintf(f, "                  int ldX,\n");
  fprintf(f, "                  double *Gi,\n");
  fprintf(f, "                  double *Gj,\n");
  fprintf(f, "                  int stG, \n");
  fprintf(f, "                  int ldG, \n");
  fprintf(f, "                  double *weights) {\n");	 

  int partial_size = 0;
  for(int i = 0; i < lA; ++i) {
    partial_size += (i + 1) * (i + 2) / 2;
  }

  fprintf(f, "   double temp[%d * NPTS_LOCAL];\n\n", size - partial_size);
  fprintf(f, "   for(int i = 0; i < %d * NPTS_LOCAL; ++i) {\n", size - partial_size);
  fprintf(f, "      temp[i] = 0.0;\n");
  fprintf(f, "   }\n\n");

  fprintf(f, "   double X_AB = shpair.rAB.x;\n");
  fprintf(f, "   double Y_AB = shpair.rAB.y;\n");
  fprintf(f, "   double Z_AB = shpair.rAB.z;\n");
  fprintf(f, "\n");

  fprintf(f, "   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {\n");
  fprintf(f, "      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);\n");
  fprintf(f, "      point *_point_outer = (_points + p_outer);\n\n");
//fprintf(f, "      double xA = shellA.origin.x;\n");
//fprintf(f, "      double yA = shellA.origin.y;\n");
//fprintf(f, "      double zA = shellA.origin.z;\n\n");
//fprintf(f, "      double xA = shpair.rA.x;\n");
//fprintf(f, "      double yA = shpair.rA.y;\n");
//fprintf(f, "      double zA = shpair.rA.z;\n\n");

//fprintf(f, "      double xB = shellB.origin.x;\n");
//fprintf(f, "      double yB = shellB.origin.y;\n");
//fprintf(f, "      double zB = shellB.origin.z;\n");
//fprintf(f, "      double xB = shpair.rB.x;\n");
//fprintf(f, "      double yB = shpair.rB.y;\n");
//fprintf(f, "      double zB = shpair.rB.z;\n");
//fprintf(f, "\n");
//fprintf(f, "      double X_AB = (xA - xB);\n");
//fprintf(f, "      double Y_AB = (yA - yB);\n");
//fprintf(f, "      double Z_AB = (zA - zB);\n");
//fprintf(f, "      double beta_in = 0.0;\n");
//fprintf(f, "      for(int i = 0; i < shellA.m; ++i) {\n");
//fprintf(f, "         for(int j = 0; j < shellB.m; ++j) {\n");
  fprintf(f, "      for(int ij = 0; ij < shpair.nprim_pair; ++ij ) {\n");
//fprintf(f, "         double aA = shellA.coeff[i].alpha;\n");
//fprintf(f, "         double cA = shellA.coeff[i].coeff;\n");
//fprintf(f, "\n");
//fprintf(f, "         double aB = shellB.coeff[j].alpha;\n");
//fprintf(f, "         double cB = shellB.coeff[j].coeff;\n");
//fprintf(f, "\n");
//fprintf(f, "         double RHO = aA + aB;\n");
  fprintf(f, "         double RHO = shpair.prim_pairs[ij].gamma;\n");
  fprintf(f, "         double RHO_INV = 1.0 / RHO;\n");
  fprintf(f, "\n");
//fprintf(f, "         double xP = (aA * xA + aB * xB) * RHO_INV;\n");
//fprintf(f, "         double yP = (aA * yA + aB * yB) * RHO_INV;\n");
//fprintf(f, "         double zP = (aA * zA + aB * zB) * RHO_INV;\n");
  fprintf(f, "         double xP = shpair.prim_pairs[ij].P.x;\n");
  fprintf(f, "         double yP = shpair.prim_pairs[ij].P.y;\n");
  fprintf(f, "         double zP = shpair.prim_pairs[ij].P.z;\n");
  fprintf(f, "\n");  
//fprintf(f, "         double X_PA = (xP - xA);\n");
//fprintf(f, "         double Y_PA = (yP - yA);\n");
//fprintf(f, "         double Z_PA = (zP - zA);\n");
  fprintf(f, "         double X_PA = shpair.prim_pairs[ij].PA.x;\n");
  fprintf(f, "         double Y_PA = shpair.prim_pairs[ij].PA.y;\n");
  fprintf(f, "         double Z_PA = shpair.prim_pairs[ij].PA.z;\n");
  fprintf(f, "\n");
//fprintf(f, "         double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);\n");
  fprintf(f, "         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;\n");
  fprintf(f, "\n");
  fprintf(f, "         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {\n");
  fprintf(f, "            point C = *(_point_outer + p_inner);\n");
  fprintf(f, "\n");
  fprintf(f, "            double xC = C.x;\n");
  fprintf(f, "            double yC = C.y;\n");
  fprintf(f, "            double zC = C.z;\n");
  fprintf(f, "\n");
  fprintf(f, "            double X_PC = (xP - xC);\n");
  fprintf(f, "            double Y_PC = (yP - yC);\n");
  fprintf(f, "            double Z_PC = (zP - zC);\n");
  fprintf(f, "\n");
  fprintf(f, "            double ");
  for(int l = 0; l < (lA + lB); ++l) {
    for(int k = 0; k < (lA + lB + 1) - l; ++k) {
      fprintf(f, "t%d%d, ", l, k);
    }
  }
  fprintf(f, "t%d%d;\n", (lA + lB), 0);
  fprintf(f, "\n");
  fprintf(f, "            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);\n");
  fprintf(f, "\n");
  traverse_dfs_vrr(f, lA, lB, root_node);
//fprintf(f, "            beta_in = 1.0;\n");
  fprintf(f, "         }\n");
  fprintf(f, "      }\n");
  fprintf(f, "\n");
  fprintf(f, "      for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {\n");
  fprintf(f, "         double *Xik = (Xi + (NPTS_LOCAL * p_outer + p_inner) * stX);\n");
  fprintf(f, "         double *Xjk = (Xj + (NPTS_LOCAL * p_outer + p_inner) * stX);\n");
  fprintf(f, "         double *Gik = (Gi + (NPTS_LOCAL * p_outer + p_inner) * stG);\n");
  fprintf(f, "         double *Gjk = (Gj + (NPTS_LOCAL * p_outer + p_inner) * stG);\n");
  fprintf(f, "\n");
  
  if(type == 0) {
    fprintf(f, "         for(int c0 = 0; c0 <= %d; ++c0) {\n", lB);
    fprintf(f, "            for(int c1 = 0; c1 <= c0; ++c1) {\n");
    fprintf(f, "               int m = %d - c0;\n", lB);
    fprintf(f, "               int n = c0 - c1;\n");
    fprintf(f, "               int p = c1;\n");
    fprintf(f, "\n");
    fprintf(f, "               int idxB = (((%d - m) * (%d - m + 1)) >> 1) + p;\n", lB, lB);
    fprintf(f, "\n");
    fprintf(f, "               double X_ABp = 1.0, comb_m_i = 1.0;\n");
    fprintf(f, "               for(int i = 0; i <= m; ++i) {\n");
    fprintf(f, "                  double rcp_i;\n");
    fprintf(f, "\n");
    fprintf(f, "                  double Y_ABp = 1.0, comb_n_j = 1.0;\n");
    fprintf(f, "                  for(int j = 0; j <= n; ++j) {\n");
    fprintf(f, "                     double rcp_j;\n");
    fprintf(f, "\n");
    fprintf(f, "                     double Z_ABp = 1.0, comb_p_k = 1.0;\n");
    fprintf(f, "                     for(int k = 0; k <= p; ++k) {\n");
    fprintf(f, "                        double rcp_k;\n");
    fprintf(f, "                        int mv, pv, Lv = %d - i - j - k;\n", lA + lB);
    fprintf(f, "\n");
    fprintf(f, "                        int offset = (Lv * (Lv + 1) * (Lv + 2) - %d) / 6;\n", lA * (lA + 1) * (lA + 2));
    fprintf(f, "                        double const_value = *(weights + NPTS_LOCAL * p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;\n");
    
    int count = 0;
    for(int r0 = 0; r0 <= lA; ++r0) {
      for(int r1 = 0; r1 <= r0; ++r1) {
	int a = lA - r0;
	int c = r1;

	int idxA = index_calculation(a, c, lA);
	fprintf(f, "                        mv = %d + m - i; pv = %d + p - k;\n", a, c);
	fprintf(f, "                        double t%d = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;\n", count);
	fprintf(f, "                        *(Gik + %d * ldG) += *(Xjk + idxB * ldX) * t%d;\n", idxA, count);
	fprintf(f, "                        *(Gjk + idxB * ldG) += *(Xik + %d * ldX) * t%d;\n", idxA, count);
      
	//if (idxA != ((lA + 1) * (lA + 2) / 2 - 1)) fprintf(f, "\n");
	count++;		
      }
    }
    fprintf(f, "\n");
    fprintf(f, "                        Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * (k + 1)); comb_p_k = (comb_p_k * (p - k)) * rcp_k;\n");
    fprintf(f, "                     }\n");
    fprintf(f, "\n");
    fprintf(f, "                     Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * (j + 1)); comb_n_j = (comb_n_j * (n - j)) * rcp_j;\n");
    fprintf(f, "                  }\n");
    fprintf(f, "\n");
    fprintf(f, "                  X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * (i + 1)); comb_m_i = (comb_m_i * (m - i)) * rcp_i;\n");
    fprintf(f, "               }\n");
    fprintf(f, "            }\n");
    fprintf(f, "         }\n");
    fprintf(f, "      }\n");
  } else if (type == 1) {
    fprintf(f, "         double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;\n");
    

    int count = 0;
    fprintf(f, "         double ");
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

	fprintf(f, "         X_ABp = 1.0; comb_m_i = 1.0;\n");
	for(int i = 0; i <= m; ++i) {
	  fprintf(f, "         Y_ABp = 1.0; comb_n_j = 1.0;\n");
	  for(int j = 0; j <= n; ++j) {
	    fprintf(f, "         Z_ABp = 1.0; comb_p_k = 1.0;\n");
	    for(int k = 0; k <= p; ++k) {
	      fprintf(f, "         const_value = *(weights + p_outer * NPTS_LOCAL + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;\n");

	      int count = 0;
	      for(int r0 = 0; r0 <= lA; ++r0) {
		for(int r1 = 0; r1 <= r0; ++r1) {
		  int a = lA - r0;
		  int c = r1;

		  int idxA = index_calculation(a, c, lA);

		  int idx = index_calculation(a + m - i, c + p - k, lA + lB - i - j - k);

		  int LAB = lA + lB - i - j - k;
		  int offset = (LAB * (LAB + 1) * (LAB + 2) - lA * (lA + 1) * (lA + 2)) / 6;
		  
		  fprintf(f, "         t%d = *(temp + %d * NPTS_LOCAL + p_inner) * const_value;\n", count, offset + idx);
		  fprintf(f, "         *(Gik + %d * ldG) += *(Xjk + %d * ldX) * t%d;\n", idxA, idxB, count);
		  fprintf(f, "         *(Gjk + %d * ldG) += *(Xik + %d * ldX) * t%d;\n", idxB, idxA, count);
      
		  count++;		
		}
	      }
	      
	      if(k < p) {
		fprintf(f, "         Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * %d); comb_p_k = (comb_p_k * %d) * rcp_k;\n", k + 1, p - k);
	      }
	    }

	    if(j < n) {
	      fprintf(f, "         Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * %d); comb_n_j = (comb_n_j * %d) * rcp_j;\n", j + 1, n - j);
	    }
	  }

	  if(i < m) {
	    fprintf(f, "         X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * %d); comb_m_i = (comb_m_i * %d) * rcp_i;\n", i + 1, m - i);
	  }
	}
	//if(idxB != ((lB + 1) * (lB + 2) / 2 - 1)) fprintf(f, "\n");
      }
    }
    fprintf(f, "      }\n");
  } else {
    fprintf(f, "Type not defined\n");
  }  
  
  fprintf(f, "   }\n");
  fprintf(f, "}\n");
}

void generate_diagonal_header_files(int lA) {
  char filename[512];
      
  sprintf(filename, "integral_%d.h", lA);
      
  FILE *f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_%d\n", lA);
  fprintf(f, "#define __MY_INTEGRAL_%d\n", lA);
  fprintf(f, "\n");
  fprintf(f, "#include \"integral_%d.h\"\n", lA);
  fprintf(f, "\n");
  fprintf(f, "void integral_%d(size_t npts,\n", lA);
//fprintf(f, "               shells shellA,\n");
  fprintf(f, "               shell_pair shpair,\n");
  fprintf(f, "               point *points,\n");
  fprintf(f, "               double *Xi,\n");
  fprintf(f, "               int stX,\n");
  fprintf(f, "               int ldX,\n");	 
  fprintf(f, "               double *Gi,\n");
  fprintf(f, "               int stG, \n");
  fprintf(f, "               int ldG, \n");
  fprintf(f, "               double *weights);\n");	   
  fprintf(f, "\n");
  fprintf(f, "#endif\n");
  
  fclose(f);
}

void generate_off_diagonal_header_files(int lA, int lB) {
  char filename[512];
      
  sprintf(filename, "integral_%d_%d.h", lA, lB);
      
  FILE *f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_%d_%d\n", lA, lB);
  fprintf(f, "#define __MY_INTEGRAL_%d_%d\n", lA, lB);
  fprintf(f, "\n");
  fprintf(f, "#include \"integral_%d_%d.h\"\n", lA, lB);
  fprintf(f, "\n");
  fprintf(f, "void integral_%d_%d(size_t npts,\n", lA, lB);
//fprintf(f, "                  shells shellA,\n");
//fprintf(f, "                  shells shellB,\n");
  fprintf(f, "                  shell_pair shpair,\n");
  fprintf(f, "                  point *points,\n");
  fprintf(f, "                  double *Xi,\n");
  fprintf(f, "                  double *Xj,\n");
  fprintf(f, "                  int stX,\n");
  fprintf(f, "                  int ldX,\n");	 
  fprintf(f, "                  double *Gi,\n");
  fprintf(f, "                  double *Gj,\n");
  fprintf(f, "                  int stG, \n");
  fprintf(f, "                  int ldG, \n");
  fprintf(f, "                  double *weights);\n");	   
  fprintf(f, "\n");
  fprintf(f, "#endif\n");
  
  fclose(f);
}

void generate_main_files(int lA) {
  char filename[512];

  FILE *f;
  
  sprintf(filename, "obara_saika_integrals.h");
      
  f = fopen(filename, "w");

  fprintf(f, "#ifndef __MY_INTEGRAL_OBARA_SAIKA\n");
  fprintf(f, "#define __MY_INTEGRAL_OBARA_SAIKA\n");
  fprintf(f, "\n");
  fprintf(f, "void compute_integral_shell_pair(size_t npts,\n");
  fprintf(f, "                  int i,\n");
  fprintf(f, "                  int j,\n");
  fprintf(f, "                  shells *shell_list,\n");
  fprintf(f, "                  point *points,\n");
  fprintf(f, "                  double *Xi,\n");
  fprintf(f, "                  double *Xj,\n");
  fprintf(f, "                  int stX,\n");
  fprintf(f, "                  int ldX,\n");	 
  fprintf(f, "                  double *Gi,\n");
  fprintf(f, "                  double *Gj,\n");
  fprintf(f, "                  int stG, \n");
  fprintf(f, "                  int ldG, \n");
  fprintf(f, "                  double *weights);\n");	   
  fprintf(f, "\n");
  fprintf(f, "#endif\n");
  
  fclose(f);  

  sprintf(filename, "obara_saika_integrals.cxx");
      
  f = fopen(filename, "w");

  fprintf(f, "#include <stdio.h>\n");
  fprintf(f, "#include <stdlib.h>\n");
  fprintf(f, "#include \"integral_data_types.h\"\n");
  fprintf(f, "#include \"obara_saika_integrals.h\"\n");
  for(int i = 0; i <= lA; ++i) {
    fprintf(f, "#include \"integral_%d.h\"\n", i);
  }

  for(int i = 0; i <= lA; ++i) {
    for(int j = 0; j <= i; ++j) {
      fprintf(f, "#include \"integral_%d_%d.h\"\n", i, j);
    }
  }
  
  fprintf(f, "\n");
  fprintf(f, "void compute_integral_shell_pair(size_t npts,\n");
  fprintf(f, "                  int i,\n");
  fprintf(f, "                  int j,\n");
  fprintf(f, "                  shells *shell_list,\n");
  fprintf(f, "                  point *points,\n");
  fprintf(f, "                  double *Xi,\n");
  fprintf(f, "                  double *Xj,\n");
  fprintf(f, "                  int stX,\n");
  fprintf(f, "                  int ldX,\n");	 
  fprintf(f, "                  double *Gi,\n");
  fprintf(f, "                  double *Gj,\n");
  fprintf(f, "                  int stG, \n");
  fprintf(f, "                  int ldG, \n");
  fprintf(f, "                  double *weights) {\n");

  fprintf(f, "   shell_pair shpair;\n");
  fprintf(f, "   // Account for permutational symmetry in kernels\n");
  fprintf(f, "   if( shell_list[i].L >= shell_list[j].L )\n");
  fprintf(f, "     generate_shell_pair(shell_list[i], shell_list[j], shpair);\n");
  fprintf(f, "   else\n");
  fprintf(f, "     generate_shell_pair(shell_list[j], shell_list[i], shpair);\n\n");
  fprintf(f, "   if (i == j) {\n");
  fprintf(f, "      int lA = shell_list[i].L;\n");
  fprintf(f, "\n");
  fprintf(f, "      if(lA == %d) {\n", 0);
  fprintf(f, "         integral_%d(npts,\n", 0);
//fprintf(f, "                    shell_list[i],\n");
  fprintf(f, "                    shpair,\n");
  fprintf(f, "                    points,\n");
  fprintf(f, "                    Xi,\n");
  fprintf(f, "                    stX,\n");
  fprintf(f, "                    ldX,\n");
  fprintf(f, "                    Gi,\n");
  fprintf(f, "                    stG, \n");
  fprintf(f, "                    ldG, \n");
  fprintf(f, "                    weights);\n");	   
  fprintf(f, "      } else ");

  for(int i = 1; i <= lA; ++i) {
    fprintf(f, "if(lA == %d) {\n", i);
    fprintf(f, "        integral_%d(npts,\n", i);
//  fprintf(f, "                   shell_list[i],\n");
    fprintf(f, "                   shpair,\n");
    fprintf(f, "                   points,\n");
    fprintf(f, "                   Xi,\n");
    fprintf(f, "                   stX,\n");
    fprintf(f, "                   ldX,\n");
    fprintf(f, "                   Gi,\n");
    fprintf(f, "                   stG, \n");
    fprintf(f, "                   ldG, \n");
    fprintf(f, "                   weights);\n");
    fprintf(f, "      } else ");
  }

  fprintf(f, "{\n");
  fprintf(f, "         printf(\"Type not defined!\\n\");\n");
  fprintf(f, "      }\n");  
  fprintf(f, "   } else {\n");
  fprintf(f, "      int lA = shell_list[i].L;\n");
  fprintf(f, "      int lB = shell_list[j].L;\n");
  fprintf(f, "\n");
  fprintf(f, "      if((lA == %d) && (lB == %d)) {\n", 0, 0);
  fprintf(f, "         integral_%d_%d(npts,\n", 0, 0);
//fprintf(f, "                      shell_list[i],\n");
//fprintf(f, "                      shell_list[j],\n");
  fprintf(f, "                      shpair,\n");
  fprintf(f, "                      points,\n");
  fprintf(f, "                      Xi,\n");
  fprintf(f, "                      Xj,\n");
  fprintf(f, "                      stX,\n");
  fprintf(f, "                      ldX,\n");
  fprintf(f, "                      Gi,\n");
  fprintf(f, "                      Gj,\n");
  fprintf(f, "                      stG, \n");
  fprintf(f, "                      ldG, \n");
  fprintf(f, "                      weights);\n");	   
  fprintf(f, "      } else ");

  for(int i = 1; i <= lA; ++i) {
    for(int j = 0; j < i; ++j) {
      fprintf(f, "if((lA == %d) && (lB == %d)) {\n", i, j);
      fprintf(f, "            integral_%d_%d(npts,\n", i, j);
    //fprintf(f, "                         shell_list[i],\n");
    //fprintf(f, "                         shell_list[j],\n");
      fprintf(f, "                         shpair,\n");
      fprintf(f, "                         points,\n");
      fprintf(f, "                         Xi,\n");
      fprintf(f, "                         Xj,\n");
      fprintf(f, "                         stX,\n");
      fprintf(f, "                         ldX,\n");
      fprintf(f, "                         Gi,\n");
      fprintf(f, "                         Gj,\n");
      fprintf(f, "                         stG, \n");
      fprintf(f, "                         ldG, \n");
      fprintf(f, "                         weights);\n");
      fprintf(f, "      } else if((lA == %d) && (lB == %d)) {\n", j, i);
      fprintf(f, "         integral_%d_%d(npts,\n", i, j);
    //fprintf(f, "                      shell_list[j],\n");
    //fprintf(f, "                      shell_list[i],\n");
      fprintf(f, "                      shpair,\n");
      fprintf(f, "                      points,\n");
      fprintf(f, "                      Xj,\n");
      fprintf(f, "                      Xi,\n");
      fprintf(f, "                      stX,\n");
      fprintf(f, "                      ldX,\n");
      fprintf(f, "                      Gj,\n");
      fprintf(f, "                      Gi,\n");
      fprintf(f, "                      stG, \n");
      fprintf(f, "                      ldG, \n");
      fprintf(f, "                      weights);\n");
      fprintf(f, "      } else ");
    }

    fprintf(f, "if((lA == %d) && (lB == %d)) {\n", i, i);
    fprintf(f, "        integral_%d_%d(npts,\n", i, i);
  //fprintf(f, "                     shell_list[i],\n");
  //fprintf(f, "                     shell_list[j],\n");
    fprintf(f, "                     shpair,\n");
    fprintf(f, "                     points,\n");
    fprintf(f, "                     Xi,\n");
    fprintf(f, "                     Xj,\n");
    fprintf(f, "                     stX,\n");
    fprintf(f, "                     ldX,\n");
    fprintf(f, "                     Gi,\n");
    fprintf(f, "                     Gj,\n");
    fprintf(f, "                     stG, \n");
    fprintf(f, "                     ldG, \n");
    fprintf(f, "                     weights);\n");
    fprintf(f, "      } else ");
  }

  fprintf(f, "{\n");
  fprintf(f, "         printf(\"Type not defined!\\n\");\n");
  fprintf(f, "      }\n");
  fprintf(f, "   }\n");  
  fprintf(f, "  delete shpair.prim_pairs;\n" );
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
      
	sprintf(filename, "integral_%d.cxx", i);
      
	FILE *f = fopen(filename, "w");

	generate_diagonal_files(f, i, size, node_list, type);

	fclose(f);
      }

      char filename[512];
      
      sprintf(filename, "integral_%d_%d.cxx", i, j);
      
      FILE *f = fopen(filename, "w");

      generate_off_diagonal_files(f, i, j, size, node_list, type);

      fclose(f);
  
      free(node_list);
    }
  }
  
  return 0;
}
