#include <vector>  // vector functionality
#include <memory>  // smart pointers
#include <fstream> // file output
#include <iostream> // terminal in and output (debugging)
#include <numeric> // for std::accumulate

#include <Readfile.h>
#include <Einsteinvisco.h>

//constructor
EinsteinVisco::EinsteinVisco(std::unique_ptr<Readfile> const &ptrReadfile)
{
  m_volume = ptrReadfile->get_volume();
  m_temperature = ptrReadfile->get_temperature();
  m_offdiagoriginal = ptrReadfile->get_transposedoffdiag();
  m_time = ptrReadfile->get_time();
  m_timestep = ptrReadfile->calc_timestep();
}


//average volume
void EinsteinVisco::m_averagevolume()
{
  m_volavg = std::accumulate(m_volume.begin(), m_volume.end(), 0.0) / m_volume.size();
}

//average temperature
void EinsteinVisco::m_averagetemperature()
{
  m_tempavg = std::accumulate(m_temperature.begin(), m_temperature.end(), 0.0) / m_temperature.size();
}

//transpose symmetrized P tensor
//void EinsteinVisco::m_transpose_symmetrized(std::vector<std::vector<double>> sym)
//{
//  for (int i=0; i<3; ++i)
//  {
//    std::vector<double> temptranspose;
//    for (std::size_t j=0; j < sym.size(); ++j)
//    { 
//      temptranspose.push_back(sym[j][i]);
//    }
//    m_transposesym.push_back(temptranspose);
//  }
//}

//symmetrize P tensor
void EinsteinVisco::m_symmetrize()
{
  std::vector<std::vector<double>> symoffdiag;
  //symmetrize tensor and store in a new array
  for (std::size_t i=0; i < m_offdiagoriginal[0].size(); ++i)
  {
    std::vector<double> offdiagvectortemp;
    //XY
    offdiagvectortemp.push_back(0.5 * (m_offdiagoriginal[0][i] + m_offdiagoriginal[2][i]));
    //XZ
    offdiagvectortemp.push_back(0.5 * (m_offdiagoriginal[1][i] + m_offdiagoriginal[4][i]));
    //YZ
    offdiagvectortemp.push_back(0.5 * (m_offdiagoriginal[3][i] + m_offdiagoriginal[5][i]));
    symoffdiag.push_back(offdiagvectortemp);
  }
  std::cout << "new stuff has number of lines: " << symoffdiag.size() << std::endl;
  
  m_transposesym = symoffdiag;
}

void EinsteinVisco::m_integrate()
{
  //integrate offdiagonal components and square the result
  for (std::size_t i=0; i<m_transposesym.size(); ++i)
  {
    std::vector<double> offdiaginttemp;
    offdiaginttemp.push_back(0.0);
    for (std::size_t j=1; j < m_transposesym[0].size(); ++j)
    {
      offdiaginttemp.push_back(offdiaginttemp[j-1] + m_transposesym[i][j] * m_timestep);
    } 
    m_alloffdiagint.push_back(offdiaginttemp);
  }
  std::cout << "Integration successful. Integration Vector size:  " << m_alloffdiagint[0].size() << std::endl;
}

void EinsteinVisco::m_transpose_integral()
{
  std::vector<std::vector<double>> transposedalloffdiagint;
  for (int i=0; i<6; ++i)
  {
    std::vector<double> temptranspose;
    for (std::size_t j=0; j < m_alloffdiagint.size(); ++j)
    {
      temptranspose.push_back(m_alloffdiagint[j][i]);
    }
    transposedalloffdiagint.push_back(temptranspose);
  }
  m_alloffdiagint = transposedalloffdiagint;
}


void EinsteinVisco::m_viscosity()
{
  size_t lines = m_alloffdiagint.size();
  size_t nsets = m_alloffdiagint[0].size();
  int range = 45000;

  std::vector<double> disq;
  for(std::size_t k=0; k <= nsets; ++k){disq.push_back(0);}

  for(std::size_t i=0; i <= lines - range; ++i)
  {
    for(std::size_t k=0; k <= nsets; ++k){disq[k] = 0;}
    std::cout << i << std::endl;
    for(std::size_t j=0; j < lines - i; ++j)
    {
      //std::cout << "j = " << j << std::endl;
      for(std::size_t k=0; k <= nsets; ++k)
      { 
        //di.push_back(m_alloffdiagint[k][j+1] - m_alloffdiagint[k][j]);
        double di = m_alloffdiagint[j+i][k] - m_alloffdiagint[j][k];
        disq[k] += di * di;
      }
    }
  }
}


//calculate einstein visco
void EinsteinVisco::calcevisco()
{
  m_averagevolume();
  m_averagetemperature();
  m_symmetrize();
  m_integrate();
  //m_transpose_integral();
  m_viscosity();
  

  std::cout << "avg. volume in nm^3   :  " << m_volavg << std::endl;
  std::cout << "avg. temperature in K  :  " << m_tempavg << std::endl;
}
