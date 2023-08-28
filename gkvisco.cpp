#include <vector>  // vector functionality
#include <memory>  // smart pointers
#include <fstream> // file output
#include <iostream> // terminal in and output (debugging)
#include <numeric> // for std::accumulate
#include <math.h>  // for square root

#include <fftw3.h>
#include <Readfile.h>
#include <GKvisco.h>

# define M_PI           3.14159265358979323846  /* pi */
# define KB             1.380649 * 1e-23 /* boltzmann constant */

//constructor
GKvisco::GKvisco(std::unique_ptr<Readfile> const &ptrReadfile)
{
  m_volume = ptrReadfile->get_volume();
  m_temperature = ptrReadfile->get_temperature();
  m_offdiagoriginal = ptrReadfile->get_transposedoffdiag();
  m_time = ptrReadfile->get_time();
  m_timestep = ptrReadfile->calc_timestep();
}


/*//average volume
void GKvisco::m_averagevolume()
{
  m_volavg = std::accumulate(m_volume.begin(), m_volume.end(), 0.0) / m_volume.size();
}
*/

void GKvisco::m_averagevolume()
{
  m_volavg = m_volume;
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
  std::vector<double> offdiagvectortemp1, offdiagvectortemp2, offdiagvectortemp3, offdiagvectortemp4, offdiagvectortemp5, offdiagvectortemp6;
  for (std::size_t i=0; i < m_offdiagoriginal[0].size(); ++i)
  {
    //XY
    offdiagvectortemp1.push_back(0.5 * (m_offdiagoriginal[0][i] + m_offdiagoriginal[2][i]));
    //XZ
    offdiagvectortemp2.push_back(0.5 * (m_offdiagoriginal[1][i] + m_offdiagoriginal[4][i]));
    //YZ
    offdiagvectortemp3.push_back(0.5 * (m_offdiagoriginal[3][i] + m_offdiagoriginal[5][i]));
    //(XX-YY)/2
    offdiagvectortemp4.push_back(0.5 * (m_offdiagoriginal[6][i] - m_offdiagoriginal[7][i]));
    //(XX-ZZ)/2
    offdiagvectortemp5.push_back(0.5 * (m_offdiagoriginal[6][i] - m_offdiagoriginal[8][i]));
    //(YY-ZZ)/2
    offdiagvectortemp6.push_back(0.5 * (m_offdiagoriginal[7][i] - m_offdiagoriginal[8][i]));
  }
  symoffdiag.push_back(offdiagvectortemp1);
  symoffdiag.push_back(offdiagvectortemp2);
  symoffdiag.push_back(offdiagvectortemp3);
  symoffdiag.push_back(offdiagvectortemp4);
  symoffdiag.push_back(offdiagvectortemp5);
  symoffdiag.push_back(offdiagvectortemp6);

  //std::cout << "new stuff has number of lines: " << symoffdiag.size() << std::endl;
  
  m_sym = symoffdiag;
}

void GKvisco::m_integrate(int N)
{
  //integrate autocorrelation
  //std::cout << "Before factor" << std::endl;
  double factor = (m_volavg * 1e-26 / (KB * m_tempavg)) * m_timestep;
  //std::cout << "After factor" << std::endl;

  std::ofstream outvisco;
  outvisco.open("SELFvisco.xvg"); 

  std::vector<double> inttemp;
  inttemp.push_back(0.0);

  double integrate1, integrate2, integrate3, integrate4, integrate5, integrate6, integrateavg = 0;
  
  for (int j=1; j < N/2; ++j)
  {
    integrate1 += 0.5 * (m_fftw_data[0][j-1] + m_fftw_data[0][j]) * factor;
    integrate2 += 0.5 * (m_fftw_data[1][j-1] + m_fftw_data[1][j]) * factor;
    integrate3 += 0.5 * (m_fftw_data[2][j-1] + m_fftw_data[2][j]) * factor;
    integrate4 += 0.5 * (m_fftw_data[3][j-1] + m_fftw_data[3][j]) * factor;
    integrate5 += 0.5 * (m_fftw_data[4][j-1] + m_fftw_data[4][j]) * factor;
    integrate6 += 0.5 * (m_fftw_data[5][j-1] + m_fftw_data[5][j]) * factor;
    integrateavg += 0.5 * (m_fftw_data[0][j-1] + m_fftw_data[0][j] + m_fftw_data[1][j-1] + m_fftw_data[1][j] + m_fftw_data[2][j-1] + m_fftw_data[2][j] + m_fftw_data[3][j-1] + m_fftw_data[3][j] + m_fftw_data[4][j-1] + m_fftw_data[4][j] + m_fftw_data[5][j-1] + m_fftw_data[5][j]) / 6 * factor;
    outvisco << j*m_timestep << "      " << integrate1 << "      " << integrate2 << "     " << integrate3 << "      " << integrate4 << "      " << integrate5 << "      " << integrate6 << "      " << integrateavg << std::endl; 
  } 
  std::cout << "Integration successful." << std::endl;
}

void GKvisco::m_do_wk(int N, fftw_complex *in)
{
    fftw_complex *out;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    //std::cout<<"fftw function call start. In is     "<< in[0][0]<<std::endl;
    //std::cout<<"fftw function call start. Out is     "<< out[0][0]<<std::endl;

    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);

    for (int i = 0; i < N; ++i)
        {
            double re = out[i][0];
            double im = out[i][1];
            in[i][0] = (re*re + im*im) / static_cast<double> (N);
            in[i][1] = 0.0;
        }

    fftw_plan p2 = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p2); /* repeat as needed */
    fftw_destroy_plan(p2);
    std::vector<double> temp;
    // double norm = out[0][0] / static_cast<double>(N);
    for (int i = 0; i < N; ++i)
        {
            const double re = out[i][0]/static_cast<double>(N); // / norm;
            temp.push_back(re);
        }
    m_fftw_data.push_back(temp);

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
  double* conv3 = m_sym[3].data();
  double* conv4 = m_sym[4].data();
  double* conv5 = m_sym[5].data();
  
  int N = m_sym[0].size();
  //std::cout<<"m_sym[0].size()      "<<N<<std::endl;

  fftw_complex *conv0c, *conv1c, *conv2c, *conv3c, *conv4c, *conv5c;
  conv0c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv1c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv3c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv4c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  conv5c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);


  for (int i = 0; i < N; ++i)
  {
    conv0c[i][0] = conv0[i];
    conv0c[i][1] = 0.0;
    conv1c[i][0] = conv1[i];
    conv1c[i][1] = 0.0;
    conv2c[i][0] = conv2[i];
    conv2c[i][1] = 0.0;
    conv3c[i][0] = conv3[i];
    conv3c[i][1] = 0.0;
    conv4c[i][0] = conv4[i];
    conv4c[i][1] = 0.0;
    conv5c[i][0] = conv5[i];
    conv5c[i][1] = 0.0;

  }

  //std::cout<<"conv0=     "<< conv0[0]<<std::endl;
  //std::cout<<"conv0c=    "<< conv0c[0][0]<<std::endl; 
  //std::cout<<"conv0cIMAG="<< conv0c[0][1]<<std::endl;

  //std::cout<<"Conversion done"<<std::endl;
  //fftw_complex *wk0, *wk1, *wk2;
  //wk0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //wk1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //wk2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  m_do_wk(N, conv0c);
  m_do_wk(N, conv1c);
  m_do_wk(N, conv2c);
  m_do_wk(N, conv3c);
  m_do_wk(N, conv4c);
  m_do_wk(N, conv5c);



  std::cout<<"all fftw done"<<std::endl;
  m_integrate(N);
  std::cout << "integration done" << std::endl;
  //m_transpose_integral();
  m_viscosity();

  //std::cout<<"0 REAL=    " <<m_fftw_data[0][0]<<std::endl;
  //std::cout<<"1 REAL=    " <<m_fftw_data[1][0]<<std::endl;
  //std::cout<<"2 REAL=    " <<m_fftw_data[2][0]<<std::endl;

  std::ofstream outac;
  outac.open("SELFautocorr.xvg");
  for (int i = 0; i<N/2; ++i)
  //for (int i = 0; i<N; ++i)
  {
    //double fftw_data_avg = (m_fftw_data[0][i] + m_fftw_data[1][i] + m_fftw_data[2][i]) / 3;
    outac << (i+1)*m_timestep << "    " << m_fftw_data[0][i] << "   " << "   " << m_fftw_data[1][i]<< "   " << m_fftw_data[2][i] << "    " << "    " << m_fftw_data[3][i] << "   " << m_fftw_data[4][i] << "   " << m_fftw_data[5][i] << std::endl;
  }

  //std::ofstream outint;
  //outint.open("SELFint.xvg");
  //for (int i = 0; i<N; ++i)
  //{
  //  outint << i*m_timestep << "      " << m_allint[0][i] << "      " << "     " << m_allint[1][i]<< "     " << m_allint[2][i] << std::endl;
  //}

  std::ofstream outvt;
  outvt.open("SELFvoltemp.xvg");
  outvt << "avg. volume in nm^3    avg. temperature in K  " << std::endl;
  outvt << m_volavg << "    " << m_tempavg << std::endl;
}
