#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <quda.h>
#include <util_quda.h>

int main()
{
  float *gauge[4];

  QudaGaugeParam Gauge_param;
  QudaInvertParam inv_param;

  Gauge_param.cpu_prec = QUDA_SINGLE_PRECISION;
  Gauge_param.cuda_prec = QUDA_SINGLE_PRECISION;
  Gauge_param.X = L1;
  Gauge_param.Y = L2;
  Gauge_param.Z = L3;
  Gauge_param.T = L4;
  Gauge_param.anisotropy = 1.0;

  Gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
  inv_param.inv_type = QUDA_CG_INVERTER;

  Gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
  Gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
  gauge_param = &Gauge_param;
  
  float mass = -0.95;
  inv_param.kappa = 1.0 / (2.0*(4 + mass));
  inv_param.tol = 1e-7;
  inv_param.maxiter = 500;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.cpu_prec = QUDA_SINGLE_PRECISION;
  inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;

  for (int dir = 0; dir < 4; dir++) {
    gauge[dir] = (float*)malloc(N*gaugeSiteSize*sizeof(float));
  }
  constructGaugeField(gauge);
  //constructUnitGaugeField(gauge);

  float *spinorIn = (float*)malloc(N*spinorSiteSize*sizeof(float));
  float *spinorOut = (float*)malloc(N*spinorSiteSize*sizeof(float));
  float *spinorCheck = (float*)malloc(N*spinorSiteSize*sizeof(float));
  if(!spinorCheck) {
    printf("malloc failed in %s\n", __func__);
    exit(1);
  }

  int i0 = 0;
  int s0 = 0;
  int c0 = 0;
  constructPointSpinorField(spinorIn, i0, s0, c0);

  double time0 = -((double)clock()); // Start the timer

  initQuda(1);
  loadQuda((void*)gauge, &Gauge_param);

  invertQuda(spinorOut, spinorIn, &inv_param);

  time0 += clock(); // stop the timer
  time0 /= CLOCKS_PER_SEC;

  printf("Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 
	 inv_param.spinorGiB, Gauge_param.gaugeGiB);
  printf("done: %i iter / %g secs = %g gflops, total time = %g secs\n", 
	 inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs, time0);

  Mat(spinorCheck, gauge, spinorOut, inv_param.kappa);
  if  (inv_param.mass_normalization == QUDA_MASS_NORMALIZATION)
    ax(0.5/inv_param.kappa, spinorCheck, N*spinorSiteSize);
  mxpy(spinorIn, spinorCheck, N*spinorSiteSize);
  float nrm2;
  nrm2 = norm(spinorCheck, N*spinorSiteSize);
  printf("r2 = %g\n", nrm2);

  endQuda();

  return 0;
}
