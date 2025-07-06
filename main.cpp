#include <iostream>
#include <vector>

#define IDX(I,J) ((I)*nx + (J))

void fillPositionVector(std::vector<double> &vec, const double& start, const double& end, const int &n)
{
  double h{((end-start)/(static_cast<double>(n)-1))};
  int i=0;
  for(auto &x :vec)
  {
    x=start+i*h;
    ++i;
  }
}

template<typename T> void printVector(const std::vector<T> &vec, const int &n)
{
  std::cout<<"Printing Vector"<<"\n";
  for(auto& x:vec)
  {
    std::cout<<x<<"\n";
  }
}

void buildBMatrix(std::vector<double> &b, const double &rho, const double& dDt, 
                  const std::vector<double> &u, const std::vector<double> &v, 
                  const double &dDx, const double &dDy, const double& nx,
                  const double& ny)
{
  double dudx;
  double dvdy;
  double dudy;
  double dvdx;
  for(int i=1;i<ny-1;++i)
  {
    for(int j=1;j<nx-1; ++j)
    {
      dudx=(u[IDX(i,j+1)] - u[IDX(i,j-1)])*0.5*dDx;
      dvdy=(v[IDX(i+1,j)] - v[IDX(i-1,j)])*0.5*dDy;
      dudy=(u[IDX(i+1,j)] - u[IDX(i-1,j)])*0.5*dDy;
      dvdx=(v[IDX(i,j+1)] - v[IDX(i,j-1)])*0.5*dDx;
      
      b[IDX(i,j)] = (rho*(dDt*(dudx+dvdy-(dudx*dudx)-2*dudy*dvdx-(dvdy*dvdy))));
    }
  }
}


void pressurePoisson(std::vector<double> &p, const double &dx, const double& dy,
                     const int &nx, const int &ny, const std::vector<double> &b,
                     const int &nit)
{
  std::vector<double> pn(p.size(),0);
  std::copy(p.begin(),p.end(),pn.begin());
  for(int q=0; q<nit;++q)
  {
    std::copy(p.begin(),p.end(),pn.begin());
    for(int i=1;i<ny-1;++i)
    {
      for(int j=1;j<nx-1;++j)
      {
        p[IDX(i,j)] = (((pn[IDX(i,j+1)] + pn[IDX(i,j-1)]*dy*dy +
                      (pn[IDX(i+1,j)] + pn[IDX(i-1,j)]*dx*dx) /
                      (2*(dx*dx + dy*dy)) - 






}


int main()
{

  // Vector definitions
  constexpr int nx{41}; // Number of points in x direction
  constexpr int ny{41}; // Number of points in y direction
  constexpr int nt{500}; // Number of time steps 
  constexpr int nit{50}; // Iteration loop for PDE?
  constexpr auto c{1}; // No clue
  
  // Time steps
  constexpr double dx{2.0/(static_cast<double>(nx)-1)}; //Difference in x vector
  constexpr double dy{2.0/(static_cast<double>(ny)-1)}; // Difference in y vector
  
  // vectors
  std::vector<double> x(nx); // x position vector
  std::vector<double> y(ny); // y position vector
  fillPositionVector(x,0,2,nx); // fill x position vector 0 to 2 with nx values
  fillPositionVector(y,0,2,ny); // fill y position vector 0 to 2 with ny values
  std::vector<double> u(nx*ny); // u vector
  std::vector<double> v(nx*ny); // v vector
  std::vector<double> p(nx*ny); // p pressure vector
  std::vector<double> b(nx*ny); // b vector?
  
  // material properties
  constexpr double rho=1.0;
  constexpr double nu=0.1;
  constexpr double dt=0.001;

  constexpr double dDx=1.0/dx;
  constexpr double dDy=1.0/dy;
  constexpr double dDt=1.0/dt;

  buildBMatrix(b,rho,dDt,u,v,dDx,dDy,nx,ny);


  return 0;
}

