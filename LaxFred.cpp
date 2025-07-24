#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>

/**
 * 1D Shock Tube 
 * C++ implementation of Lax-Friedrich and Adjustable Time Stepping
 * https://uk.mathworks.com/matlabcentral/fileexchange/82240-1d-shock-tube
**/


void initialiseX(std::vector<double> &vec,const int &N, const double &start, const double &end)
{
	double h=(end-start)/(static_cast<double>(N)-1);
	int i=0;
	for(auto &x: vec)
	{
		x=start+i*h;
		++i;
	}
}

void initialiseXC(const std::vector<double> &x, std::vector<double> &xc,const int &N)
{
	for(int i=0;i<N;++i)
	{
		xc[i]=0.5*(x[i]+x[i+1]);
	}
}

void initialiseFlow(std::vector<double> &rho, std::vector<double> &p, std::vector<double> &u,
										std::vector<double> &e, const double &rhoL, const double &pL, 
										const double& rhoR, const double & pR, const double& gamma, const int &N)
{
	for(int i=0; i<N; ++i)
	{
		if(i<=N/2)
		{
			rho[i]=rhoL;
			p[i]=pL;
		}
		else
		{
			rho[i]=rhoR;
			p[i]=pR;
		}
		u[i]=0;
	}
	for(int i=0;i<N;++i)
	{
		e[i]=(p[i]/(gamma-1))+0.5*rho[i]*u[i]*u[i];
	}
}


void LaxFriedrichs(double &time, const double& endTime, const int &N, const double& gamma, const double& CFL,
									 std::vector<double> &p, std::vector<double> &e, std::vector<double> &rho, std::vector<double> &u,
									 std::vector<double> &xc, std::vector<double> &NewRho, std::vector<double> &NewU, std::vector<double> &NewE)
{
	while(time<=endTime)
	{
		for(int i=1;i<N-1;++i)
		{
			p[i]=(gamma-1)*(e[i]-0.5*rho[i]*u[i]*u[i]);
		}
		std::vector<double> a(N);
		for(int i=0;i<N;++i)
		{
			a[i]=sqrt(gamma*p[i]/rho[i]);
		}
		double lambda=*std::max_element(a.begin(),a.end());
		double max_vel=*std::max_element(u.begin(),u.end());

		double dt=CFL/static_cast<double>(N)/(max_vel+lambda);
		double dx;
		
	
		time=time+dt;

		for(int i=1;i<N-1;++i)
		{
			dx=xc[i]-xc[i-1];

			//Momentum, Right
			double mom_R=rho[i+1]*u[i+1];
			double rho_R=rho[i+1];
			double u_R=u[i+1];
			double p_R=p[i+1];

			// Momentum, P?
			double mom_P=rho[i]*u[i];
			double rho_P=rho[i];
			double u_P=u[i];
			double p_P=p[i];

			//Momentum, L
			double mom_L=rho[i-1]*u[i-1];
			double rho_L=rho[i-1];
			double u_L=u[i-1];
			double p_L=p[i-1];


			//Velocity Flux, R, P and L
			double vel_flux_R=rho_R*u_R*u_R+p_R;
			double e_R=e[i+1];

			double vel_flux_P=rho_P*u_P*u_P+p_P;
			double e_P=e[i];

			double vel_flux_L=rho_L*u_L*u_L+p_L;
			double e_L=e[i-1];

			// Energy Flux, Density Flux, Velocity Flux states
			double energy_flux_R=u_R*(e_R+p_R);
      double energy_flux_P=u_P*(e_P+p_P);
      double energy_flux_L=u_L*(e_L+p_L);
        
      double rho_fluxR=0.5*(mom_P+mom_R)-0.5*lambda*(rho_R-rho_P);
      double rho_fluxL=0.5*(mom_P+mom_L)-0.5*lambda*(rho_P-rho_L);
        
      double vel_fluxR=0.5*(vel_flux_P+vel_flux_R)-0.5*lambda*(mom_R-mom_P) ;
      double vel_fluxL=0.5*(vel_flux_P+vel_flux_L)-0.5*lambda*(mom_P-mom_L) ;
        
      double energy_fluxR=0.5*(energy_flux_P+energy_flux_R)-0.5*lambda*(e_R-e_P);
      double energy_fluxL=0.5*(energy_flux_P+energy_flux_L)-0.5*lambda*(e_P-e_L);

			NewRho[i]=rho_P-(dt/dx)*(rho_fluxR-rho_fluxL);
			double vel_flux=mom_P-(dt/dx)*(vel_fluxR-vel_fluxL);

			NewU[i]=vel_flux/NewRho[i];
			NewE[i]=e_P-(dt/dx)*(energy_fluxR-energy_fluxL);
		}
		rho=NewRho;
		u=NewU;
		e=NewE;
		}
}


void writeSolution(const std::vector<double> &vec, std::string filename)
{	
	std::ofstream file(filename);
	for(auto x: vec)
	{
		file<<x<<'\n';
	}
}


int main()
{
	// Initialisation of Parameters
	constexpr int N=800;
	constexpr double gamma=1.4;
	constexpr double endTime=0.2;
	constexpr double CFL=0.3;

	//Grid Initialisation 
	std::vector<double> x(N);
	initialiseX(x,N,0,1);
	std::vector<double> xc(N);
	initialiseXC(x,xc,N);
	xc[0]=0;
	xc[N-1]=1;
	double t=0;

	//Initial Conditions
	constexpr double rhoR=0.125;
	constexpr double rhoL=1.0;
	constexpr double pR=0.1;
	constexpr double pL=1.0;
	std::vector<double> rho(N);
	std::vector<double> p(N);
	std::vector<double> u(N);
	std::vector<double> e(N);
	std::vector<double> NewRho(N);
	std::vector<double> NewU(N);
	std::vector<double> NewE(N);
	initialiseFlow(rho,p,u,e,rhoL,pL,rhoR,pR,gamma,N);

	//Copy Vectors
	NewRho=rho;
	NewU=u;
	NewE=e;

	//Run Simulation
	LaxFriedrichs(t,endTime,N,gamma,CFL,p,e,rho,u,xc,NewRho,NewU,NewE);
	writeSolution(xc,"xc.txt");
	writeSolution(rho,"rho.txt");
	writeSolution(u,"u.txt");
	writeSolution(e,"e.txt");


	std::cout<<"Completed Simulation!" <<std::endl;
	return 0;
}

