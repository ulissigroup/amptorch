#include <math.h>
#include <stdio.h>
#include "calculate_gmp.h"

// extern "C" int calculate_atomistic_mcsh(double** cell, double** cart, double** scale,
//                                         int* atom_i, int natoms, int* cal_atoms, int cal_num,
//                                         int** params_i, double** params_d, int nsyms,
//                                         double** symf, double** dsymf, double** dsymfa) {
    // cell: cell info of structure
    // cart: cartesian coordinates of atoms
    // scale: fractional coordinates of atoms
    // atom_i: atom type index (start with 1)
    // params_i: integer parameter for symmetry function
    //           [symmetry function type, 1st neighbor atom type, 2nd neighbor atom type]
    // params_d: double parameter for symmetry function
    //           [cutoff, param1, param2, param3]
    // natoms: # of atoms
    // nsyms: # of symmetry functions
    // symf: symmetry function vector ([# of atoms, # of symfuncs])
    // dsymf: derivative of symmetry function vector
    // originally, dsymf is 4D array (dimension: [# of atoms, # of symfuncs, # of atoms, 3])
    // in this function, we use 2D array ([# of atoms, # of symfuncs * # of atoms * 3])
    // Similarly, dsymfa is 4D array (dimension: [# of atoms, # of symfuncs, 3, 6])
    // in this function, we use 2D array ([# of atoms, # of symfuncs * 3 * 6])


    // params_d: sigma, weight, A, beta, cutoff
    // atom_gaussian: 2D array (dimension: [# atom types, # gaussian * 2]), i*2: B, i*2+1: alpha
    // ngaussian: number of gaussian for each atom type [n gaussian for type 1, n gaussian for type 2 ...]
extern "C" int calculate_gmp(double** cell, double** cart, double** scale, int* pbc_bools,
                                        int* atom_i, int natoms, int* cal_atoms, int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        // same atom
                        // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                        //     continue;

                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];

                        // if (tmp < cutoff) {
                        if (tmp_r2 < cutoff_sqr) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*4 + a] = total_shift[a];
                            nei_list_d[nneigh*4 + 3] = tmp_r2;
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_type = get_mcsh_type(params_i[m][0], params_i[m][1]);
            GMPFunction mcsh_function = get_mcsh_function(params_i[m][0], params_i[m][1]);

            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double M = 0;
            if (mcsh_type == 1){
                double m_desc[1], deriv[3];

                for (int j = 0; j < nneigh; ++j) {
                    double dMdx = 0.0, dMdy = 0.0, dMdz = 0.0;

                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc, deriv);
                        M += m_desc[0];
                        dMdx += deriv[0];
                        dMdy += deriv[1];
                        dMdz += deriv[2];
                    }
                    dMdx = dMdx * weight;
                    dMdy = dMdy * weight;
                    dMdz = dMdz * weight;

                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                    dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                    dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                    dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                }
                M = M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 2){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                }

                double miu[3], deriv[9];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                    }
                }
                M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                // M = sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3;
                double dMdx, dMdy, dMdz;
                if (fabs(M) <= 1e-8) {
                    M = 0;
                }
                else {
                    for (int j = 0; j < nneigh; ++j) {
                        dMdx = (1.0/M) * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * weight;
                        dMdy = (1.0/M) * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * weight;
                        dMdz = (1.0/M) * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * weight;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                        dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                        dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                        dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                    }
                }
                M = M * weight;
                mcsh[ii][m] += M;

                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
            }

            if (mcsh_type == 3){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu4_dxj = new double[nneigh];
                double* sum_dmiu5_dxj = new double[nneigh];
                double* sum_dmiu6_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu4_dyj = new double[nneigh];
                double* sum_dmiu5_dyj = new double[nneigh];
                double* sum_dmiu6_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                double* sum_dmiu4_dzj = new double[nneigh];
                double* sum_dmiu5_dzj = new double[nneigh];
                double* sum_dmiu6_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu4_dxj[j] = 0.0;
                    sum_dmiu5_dxj[j] = 0.0;
                    sum_dmiu6_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu4_dyj[j] = 0.0;
                    sum_dmiu5_dyj[j] = 0.0;
                    sum_dmiu6_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                    sum_dmiu4_dzj[j] = 0.0;
                    sum_dmiu5_dzj[j] = 0.0;
                    sum_dmiu6_dzj[j] = 0.0;
                }

                double miu[6], deriv[18];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_miu4 += miu[3];
                        sum_miu5 += miu[4];
                        sum_miu6 += miu[5];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                        sum_dmiu4_dxj[j] += deriv[9];
                        sum_dmiu4_dyj[j] += deriv[10];
                        sum_dmiu4_dzj[j] += deriv[11];
                        sum_dmiu5_dxj[j] += deriv[12];
                        sum_dmiu5_dyj[j] += deriv[13];
                        sum_dmiu5_dzj[j] += deriv[14];
                        sum_dmiu6_dxj[j] += deriv[15];
                        sum_dmiu6_dyj[j] += deriv[16];
                        sum_dmiu6_dzj[j] += deriv[17];
                    }
                }
                M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                         sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                double dMdx, dMdy, dMdz;
                if (fabs(M) <= 1e-8) {
                    M = 0;
                }
                else {
                    for (int j = 0; j < nneigh; ++j) {
                        dMdx = (1.0/M) * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                        sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                        sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * weight;

                        dMdy = (1.0/M) * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                        sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                        sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * weight;

                        dMdz = (1.0/M) * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                        sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                        sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * weight;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                        dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                        dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                        dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                    }
                }

                M = M * weight;
                mcsh[ii][m] += M;

                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu4_dxj;
                delete [] sum_dmiu5_dxj;
                delete [] sum_dmiu6_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu4_dyj;
                delete [] sum_dmiu5_dyj;
                delete [] sum_dmiu6_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
                delete [] sum_dmiu4_dzj;
                delete [] sum_dmiu5_dzj;
                delete [] sum_dmiu6_dzj;
            }
        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}




extern "C" int calculate_gmp_noderiv(double** cell, double** cart, double** scale, int* pbc_bools,
                                        int* atom_i, int natoms, int* cal_atoms, int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        // same atom
                        // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                        //     continue;

                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];

                        // if (tmp < cutoff) {
                        if (tmp_r2 < cutoff_sqr) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*4 + a] = total_shift[a];
                            nei_list_d[nneigh*4 + 3] = tmp_r2;
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_type = get_mcsh_type(params_i[m][0], params_i[m][1]);
            GMPFunctionNoderiv mcsh_function = get_mcsh_function_noderiv(params_i[m][0], params_i[m][1]);

            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double M = 0.0;
            if (mcsh_type == 1){
                double m_desc[1];

                for (int j = 0; j < nneigh; ++j) {

                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc);
                        M += m_desc[0];
                    }
                }
                M = M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 2){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                double miu[3];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                    }
                }
                M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                M = M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 3){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                double miu[6];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_miu4 += miu[3];
                        sum_miu5 += miu[4];
                        sum_miu6 += miu[5];
                    }
                }
                M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                         sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                M = M * weight;
                mcsh[ii][m] += M;
            }
        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}








extern "C" int calculate_gmp_square(double** cell, double** cart, double** scale, int* pbc_bools,
                                        int* atom_i, int natoms, int* cal_atoms, int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        // same atom
                        // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                        //     continue;

                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];

                        // if (tmp < cutoff) {
                        if (tmp_r2 < cutoff_sqr) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*4 + a] = total_shift[a];
                            nei_list_d[nneigh*4 + 3] = tmp_r2;
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_type = get_mcsh_type(params_i[m][0], params_i[m][1]);
            GMPFunction mcsh_function = get_mcsh_function(params_i[m][0], params_i[m][1]);

            // params_d: sigma, weight, A, beta, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double M = 0;
            if (mcsh_type == 1){
                double m_desc[1], deriv[3];

                for (int j = 0; j < nneigh; ++j) {
                    double dMdx = 0.0, dMdy = 0.0, dMdz = 0.0;

                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc, deriv);
                        M += m_desc[0];
                        dMdx += deriv[0];
                        dMdy += deriv[1];
                        dMdz += deriv[2];
                    }
                    // dMdx = dMdx * weight;
                    // dMdy = dMdy * weight;
                    // dMdz = dMdz * weight;
                    dMdx = 2.0 * M * dMdx * weight;
                    dMdy = 2.0 * M * dMdy * weight;
                    dMdz = 2.0 * M * dMdz * weight;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                    dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                    dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                    dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                }
                // M = M * weight;
                M = M * M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 2){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                }

                double miu[3], deriv[9];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                    }
                }
                // M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                M = sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3;
                double dMdx, dMdy, dMdz;
                for (int j = 0; j < nneigh; ++j) {
                    // dMdx = (1.0/M) * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * weight;
                    // dMdy = (1.0/M) * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * weight;
                    // dMdz = (1.0/M) * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * weight;
                    dMdx = 2.0 * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * weight;
                    dMdy = 2.0 * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * weight;
                    dMdz = 2.0 * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * weight;

                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                    dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                    dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                    dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                }
                M = M * weight;
                mcsh[ii][m] += M;

                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
            }

            if (mcsh_type == 3){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                double* sum_dmiu1_dxj = new double[nneigh];
                double* sum_dmiu2_dxj = new double[nneigh];
                double* sum_dmiu3_dxj = new double[nneigh];
                double* sum_dmiu4_dxj = new double[nneigh];
                double* sum_dmiu5_dxj = new double[nneigh];
                double* sum_dmiu6_dxj = new double[nneigh];
                double* sum_dmiu1_dyj = new double[nneigh];
                double* sum_dmiu2_dyj = new double[nneigh];
                double* sum_dmiu3_dyj = new double[nneigh];
                double* sum_dmiu4_dyj = new double[nneigh];
                double* sum_dmiu5_dyj = new double[nneigh];
                double* sum_dmiu6_dyj = new double[nneigh];
                double* sum_dmiu1_dzj = new double[nneigh];
                double* sum_dmiu2_dzj = new double[nneigh];
                double* sum_dmiu3_dzj = new double[nneigh];
                double* sum_dmiu4_dzj = new double[nneigh];
                double* sum_dmiu5_dzj = new double[nneigh];
                double* sum_dmiu6_dzj = new double[nneigh];
                for (int j=0; j<nneigh; j++) {
                    sum_dmiu1_dxj[j] = 0.0;
                    sum_dmiu2_dxj[j] = 0.0;
                    sum_dmiu3_dxj[j] = 0.0;
                    sum_dmiu4_dxj[j] = 0.0;
                    sum_dmiu5_dxj[j] = 0.0;
                    sum_dmiu6_dxj[j] = 0.0;
                    sum_dmiu1_dyj[j] = 0.0;
                    sum_dmiu2_dyj[j] = 0.0;
                    sum_dmiu3_dyj[j] = 0.0;
                    sum_dmiu4_dyj[j] = 0.0;
                    sum_dmiu5_dyj[j] = 0.0;
                    sum_dmiu6_dyj[j] = 0.0;
                    sum_dmiu1_dzj[j] = 0.0;
                    sum_dmiu2_dzj[j] = 0.0;
                    sum_dmiu3_dzj[j] = 0.0;
                    sum_dmiu4_dzj[j] = 0.0;
                    sum_dmiu5_dzj[j] = 0.0;
                    sum_dmiu6_dzj[j] = 0.0;
                }

                double miu[6], deriv[18];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_miu4 += miu[3];
                        sum_miu5 += miu[4];
                        sum_miu6 += miu[5];
                        sum_dmiu1_dxj[j] += deriv[0];
                        sum_dmiu1_dyj[j] += deriv[1];
                        sum_dmiu1_dzj[j] += deriv[2];
                        sum_dmiu2_dxj[j] += deriv[3];
                        sum_dmiu2_dyj[j] += deriv[4];
                        sum_dmiu2_dzj[j] += deriv[5];
                        sum_dmiu3_dxj[j] += deriv[6];
                        sum_dmiu3_dyj[j] += deriv[7];
                        sum_dmiu3_dzj[j] += deriv[8];
                        sum_dmiu4_dxj[j] += deriv[9];
                        sum_dmiu4_dyj[j] += deriv[10];
                        sum_dmiu4_dzj[j] += deriv[11];
                        sum_dmiu5_dxj[j] += deriv[12];
                        sum_dmiu5_dyj[j] += deriv[13];
                        sum_dmiu5_dzj[j] += deriv[14];
                        sum_dmiu6_dxj[j] += deriv[15];
                        sum_dmiu6_dyj[j] += deriv[16];
                        sum_dmiu6_dzj[j] += deriv[17];
                    }
                }
                // M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                //          sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                M = sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                    sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6;
                double dMdx, dMdy, dMdz;
                for (int j = 0; j < nneigh; ++j) {
                    // dMdx = (1.0/M) * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                    //                 sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                    //                 sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * weight;

                    // dMdy = (1.0/M) * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                    //                 sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                    //                 sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * weight;

                    // dMdz = (1.0/M) * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                    //                 sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                    //                 sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * weight;

                    dMdx = (2.0) * (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                    sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                    sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * weight;

                    dMdy = (2.0) * (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                    sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                    sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * weight;

                    dMdz = (2.0) * (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                    sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                    sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * weight;

                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dMdx;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dMdy;
                    dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dMdz;

                    dmcsh[ii*nmcsh + m][i*3]     -= dMdx;
                    dmcsh[ii*nmcsh + m][i*3 + 1] -= dMdy;
                    dmcsh[ii*nmcsh + m][i*3 + 2] -= dMdz;
                }
                M = M * weight;
                mcsh[ii][m] += M;

                delete [] sum_dmiu1_dxj;
                delete [] sum_dmiu2_dxj;
                delete [] sum_dmiu3_dxj;
                delete [] sum_dmiu4_dxj;
                delete [] sum_dmiu5_dxj;
                delete [] sum_dmiu6_dxj;
                delete [] sum_dmiu1_dyj;
                delete [] sum_dmiu2_dyj;
                delete [] sum_dmiu3_dyj;
                delete [] sum_dmiu4_dyj;
                delete [] sum_dmiu5_dyj;
                delete [] sum_dmiu6_dyj;
                delete [] sum_dmiu1_dzj;
                delete [] sum_dmiu2_dzj;
                delete [] sum_dmiu3_dzj;
                delete [] sum_dmiu4_dzj;
                delete [] sum_dmiu5_dzj;
                delete [] sum_dmiu6_dzj;
            }
        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}




extern "C" int calculate_gmp_square_noderiv(double** cell, double** cart, double** scale, int* pbc_bools,
                                        int* atom_i, int natoms, int* cal_atoms, int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        // same atom
                        // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                        //     continue;

                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];

                        // if (tmp < cutoff) {
                        if (tmp_r2 < cutoff_sqr) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*4 + a] = total_shift[a];
                            nei_list_d[nneigh*4 + 3] = tmp_r2;
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_type = get_mcsh_type(params_i[m][0], params_i[m][1]);
            GMPFunctionNoderiv mcsh_function = get_mcsh_function_noderiv(params_i[m][0], params_i[m][1]);

            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double M = 0.0;
            if (mcsh_type == 1){
                double m_desc[1];

                for (int j = 0; j < nneigh; ++j) {

                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc);
                        M += m_desc[0];
                    }
                }
                // M = M * weight;
                M = M * M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 2){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                double miu[3];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                    }
                }
                // M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                M = sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3;
                M = M * weight;
                mcsh[ii][m] += M;
            }

            if (mcsh_type == 3){
                double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                double miu[6];
                for (int j = 0; j < nneigh; ++j) {
                    int neigh_atom_element_index = nei_list_i[j*2];
                    int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                    double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                    for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                        double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                        mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                        // miu: miu_1, miu_2, miu_3
                        // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                        sum_miu1 += miu[0];
                        sum_miu2 += miu[1];
                        sum_miu3 += miu[2];
                        sum_miu4 += miu[3];
                        sum_miu5 += miu[4];
                        sum_miu6 += miu[5];
                    }
                }
                // M = sqrt(sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                //          sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                M = sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                    sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6;
                M = M * weight;
                mcsh[ii][m] += M;
            }
        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}

void PyInit_libmcsh(void) { } // for windows
