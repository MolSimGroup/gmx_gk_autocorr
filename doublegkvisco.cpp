#include <vector>  // vector functionality
#include <memory>  // smart pointers
#include <fstream> // file output
#include <iostream> // terminal in and output (debugging)
#include <numeric> // for std::accumulate

#include <fftw3.h>
#include <Readfile.h>
#include <GKvisco.h>

//constructor
GKvisco::GKvisco(std::unique_ptr<Readfile> const &ptrReadfile)
{
  m_volume = ptrReadfile->get_volume();
  m_temperature = ptrReadfile->get_temperature();
  m_offdiagoriginal = ptrReadfile->get_transposedoffdiag();
  m_time = ptrReadfile->get_time();
  m_timestep = ptrReadfile->calc_timestep();
}


//average volume
void GKvisco::m_averagevolume()
{
  m_volavg = std::accumulate(m_volume.begin(), m_volume.end(), 0.0) / m_volume.size();
}

//average temperature
void GKvisco::m_averagetemperature()
{
  m_tempavg = std::accumulate(m_temperature.begin(), m_temperature.end(), 0.0) / m_temperature.size();
}

//symmetrize P tensor
void GKvisco::m_symmetrize()
{
  std::vector<std::vector<double>> symoffdiag;
  //symmetrize tensor and store in a new array
  std::vector<double> offdiagvectortemp1, offdiagvectortemp2, offdiagvectortemp3;
  for (std::size_t i=0; i < m_offdiagoriginal[0].size(); ++i)
  {
    //XY
    offdiagvectortemp1.push_back(0.5 * (m_offdiagoriginal[0][i] + m_offdiagoriginal[2][i]));
    //XZ
    offdiagvectortemp2.push_back(0.5 * (m_offdiagoriginal[1][i] + m_offdiagoriginal[4][i]));
    //YZ
    offdiagvectortemp3.push_back(0.5 * (m_offdiagoriginal[3][i] + m_offdiagoriginal[5][i]));
  }
  symoffdiag.push_back(offdiagvectortemp1);
  symoffdiag.push_back(offdiagvectortemp2);
  symoffdiag.push_back(offdiagvectortemp3);
  std::cout << "new stuff has number of lines: " << symoffdiag.size() << std::endl;
  
  m_sym = symoffdiag;
}

void GKvisco::m_integrate(int N, fftw_complex *res)
{
  //integrate autocorrelation
  for (int i=0; i < N; ++i)
  {
    std::vector<double> inttemp;
    inttemp.push_back(0.0);
    for (int j=1; j < N; ++j)
    {
      inttemp.push_back(inttemp[j-1] + res[i][0] * m_timestep);
    } 
    m_allint.push_back(inttemp);
  }
  std::cout << "Integration successful. Integration Vector size:  " << m_allint[0].size() << std::endl;
}

void GKvisco::m_do_fftw(int N, fftw_complex *in, fftw_complex *out)
{
    std::cout<<"fftw function call start. In is     "<< in[0][0]<<std::endl;
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);
    std::cout<<"fftw out     " << out[0][0]<<std::endl;
}

void GKvisco::m_do_compl_multpl(int N, fftw_complex *in, fftw_complex *out)
{
    std::cout<<"N = "<<N<<std::endl;
    std::cout<<"in[0][0] = "<<in[0][0]<<std::endl;
    std::cout<<"in[0][1] = "<<in[0][1]<<std::endl;
    for (int i = 0; i < N; ++i)
        {
            out[i][0] = in[i][0] * in[i][0] + in[i][1] * in[i][1];
            out[i][1] = 0;
        }
    std::cout<<"out in multpl    "<< out[0][0]<<std::endl;
}

double* GKvisco::m_do_ifftw(int N, fftw_complex *in)
{
    double* out = (double*) fftw_malloc(sizeof(double)*N);
    fftw_plan p = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);
    return out;
}

void GKvisco::m_viscosity()
{
  //size_t lines = m_alloffdiagint.size();
  //size_t nsets = m_alloffdiagint[0].size();
}

//calculate green kubo viscosity
void GKvisco::calcvisco()
{
  m_averagevolume();
  m_averagetemperature();
  m_symmetrize();

  std::cout<<"symmetrize ok"<<std::endl;
  std::cout<<"check data"<<std::endl;

  /* convert data from std::vector to double array */
  double* conv0 = m_sym[0].data();
  double* conv1 = m_sym[1].data();
  double* conv2 = m_sym[2].data();
  
  int N = m_sym[0].size();
  std::cout<<"m_sym.size()      "<<N<<std::endl;

  fftw_complex *conv0c, *conv1c, *conv2c;
  conv0c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv1c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for (int i = 0; i < N; ++i)
  {
    conv0c[i][0] = conv0[i];
    conv0c[i][1] = 0;
    conv1c[i][0] = conv1[i];
    conv1c[i][1] = 0;
    conv2c[i][0] = conv2[i];
    conv2c[i][1] = 0;
  }

  std::cout<<"conv0=     "<< conv0[20]<<std::endl;
  std::cout<<"conv0c=    "<< conv0c[20][0]<<std::endl; 
  std::cout<<"conv0cIMAG="<< conv0c[20][1]<<std::endl; 

  std::cout<<"Conversion done"<<std::endl;
  fftw_complex *fftwout0, *fftwout1, *fftwout2;
  fftwout0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftwout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftwout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  m_do_fftw(N, conv0c, fftwout0);
  m_do_fftw(N, conv1c, fftwout1);
  m_do_fftw(N, conv2c, fftwout2);

  /* multiply fftw results with its own complex conjugate */
  fftw_complex *multpl0, *multpl1, *multpl2;
  multpl0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  multpl1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  multpl2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  m_do_compl_multpl(N, fftwout0, multpl0);
  m_do_compl_multpl(N, fftwout1, multpl1);
  m_do_compl_multpl(N, fftwout2, multpl2);

  std::cout<<"multipl0[0][0] =     "<<multpl0[0][0]<<std::endl;

  double *res0, *res1, *res2;

  res0 = (double*) fftw_malloc(sizeof(double) * N);
  res1 = (double*) fftw_malloc(sizeof(double) * N);
  res2 = (double*) fftw_malloc(sizeof(double) * N);

  res0 = m_do_ifftw(N, multpl0); 
  res1 = m_do_ifftw(N, multpl1); 
  res2 = m_do_ifftw(N, multpl2); 

  std::cout<<"all fftw done"<<std::endl;
  //m_integrate(N);
  //m_transpose_integral();
  m_viscosity();

  for (int i = 0; i<10; ++i)
  {
    std::cout<<"res0 REAL=    " <<res0[i]<<std::endl;
    std::cout<<"res1 REAL=    " <<res1[i]<<std::endl;
    std::cout<<"res2 REAL=    " <<res2[i]<<std::endl;
  }
  

  std::cout << "avg. volume in nm^3   :  " << m_volavg << std::endl;
  std::cout << "avg. temperature in K  :  " << m_tempavg << std::endl;
}
