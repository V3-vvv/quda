#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include <quda_internal.h>
#include <eigensolve_quda.h>
#include <qio_field.h>
#include <color_spinor_field.h>
#include <blas_quda.h>
#include <util_quda.h>
#include <tune_quda.h>
#include <eigen_helper.h>

#include <random_quda.h>
#include <random>

namespace quda
{
    TestSolver::TestSolver(const DiracMatrix &mat, QudaEigParam *eig_param) : EigenSolver(mat, eig_param)
    {
        getProfile().TPSTART(QUDA_PROFILE_INIT);

        // Tridiagonal/Arrow matrix
        alpha.resize(n_kr, 0.0);
        beta.resize(n_kr, 0.0);

        // Test solver specific checks
        if (eig_param->spectrum != QUDA_SPECTRUM_LR_EIG && eig_param->spectrum != QUDA_SPECTRUM_SR_EIG) {
            errorQuda("Only real spectrum type (LR or SR) can be passed to the test solver");
        }

        getProfile().TPSTOP(QUDA_PROFILE_INIT);
    }

    void TestSolver::operator()(std::vector<ColorSpinorField> &kSpace, std::vector<Complex> &evals)
    {
        // initial v0
        // RNG rng(kSpace[0], 1234);
        // auto norm = blas::norm2(kSpace[0]);
        // if (!std::isfinite(norm) || norm == 0.0) {
        //     spinorNoise(kSpace[0], rng, QUDA_NOISE_UNIFORM);
        // }

        MatrixXd V = MatrixXd::Zero(eig_param->n_kr, kSpace[0].Volume());
        VectorXd u = VectorXd::Zero(kSpace[0].Volume());

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        for(int i = 0; i < kSpace[0].Volume(); i++)
        {
            u(i) = dis(gen);
        }

        // normalization
        double norm_u = u.norm();
        V.row(0) = u / norm_u;
        
        lanczos(kSpace, 0);

        int iter = eig_param->n_kr;

        VectorXd beta_k;
        VectorXd evals_nev;
        MatrixXd evecs_nev;
        
        ritz(beta_k, evals_nev, evecs_nev);

        MatrixXd x;
        x = V * evecs_nev;

        Complex* evals_ptr = reinterpret_cast<Complex*>(evals_nev.data());

        mat(evals_ptr, V[0]);
        

    }

    void TestSolver::lanczos(std::vector<ColorSpinorField> &v, int j)
    {
        // Compute r = A * v_j - b_{j-i} * v_{j-1}
        // r = A * v_j

        chebyOp(r[0], v[j]);

        // a_j = v_j^dag * r
        alpha[j] = blas::reDotProduct(v[j], r[0]);

        // r = r - a_j * v_j
        blas::axpy(-alpha[j], v[j], r[0]);

        int start = (j > num_keep) ? j - 1 : 0;

        if (j - start > 0) {
        std::vector<double> beta_ = {beta.begin() + start, beta.begin() + j};
        for (auto & bi : beta_) bi = -bi;

        // r = r - b_{j-1} * v_{j-1}
        blas::block::axpy(beta_, {v.begin() + start, v.begin() + j}, r[0]);
        }

        // Orthogonalise r against the Krylov space
        for (int k = 0; k < 1; k++) blockOrthogonalizeHMGS(v, r, ortho_block_size, j + 1);

        // b_j = ||r||
        beta[j] = sqrt(blas::norm2(r[0]));

        // Prepare next step.
        // v_{j+1} = r / b_j
        blas::axy(1.0 / beta[j], r[0], v[j + 1]);
    }

    void TestSolver::ritz(VectorXd& beta_k, VectorXd& evals_nev, MatrixXd& evecs_nev){
        int n = alpha.size();
        MatrixXd t = MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i) 
        {
            t(i, i) = alpha[i];
            t(i, i + 1) = beta[i];
            t(i + 1, i) = beta[i];
        }

        if(beta_k.size() > 0)
        {
            t.block(n_ev, 0, 1, n_ev) = beta_k.transpose();
            t.block(0, n_ev, n_ev, 1) = beta_k;
        }

        SelfAdjointEigenSolver<MatrixXd> eigensolver(t);
        VectorXd evals = eigensolver.eigenvalues();
        MatrixXd evecs = eigensolver.eigenvectors();

        std::vector<std::pair<double, int>> indexed_evals;
        for(int i = 0; i < evals.size(); i++)
        {
            indexed_evals.push_back(std::make_pair(evals(i), i));
        }
        std::sort(indexed_evals.begin(), indexed_evals.end(), compareIndices);

        // 选择前 k 个特征值和特征向量
        Eigen::VectorXd wk(n_ev);
        Eigen::MatrixXd sk(evecs.rows(), n_ev);
        for (int i = 0; i < n_ev; ++i) {
            wk[i] = indexed_evals[i].first;
            sk.col(i) = evecs.col(indexed_evals[i].second);
        }

        evals_nev = wk;
        evecs_nev = sk;
    }

    bool TestSolver::compareIndices(const std::pair<double, int>& a, const std::pair<double, int>& b) 
    {
        return a.first < b.first;
    }
}