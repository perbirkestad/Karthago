/* 3D Convection-Diffusion equation solved for a cube with side L
Performance 50E6 on Intel(R) Xeon(R) CPU E3-1240 v3 @ 3.40GHz single core
100x100x100 cells, 4000 timesteps (82 s)
Periodic BC:s   										*/


#include <iostream>
#include <fstream>

using namespace std;

int const N=64;
int number_of_timesteps=4000;
float c[N*N*N];
float dcdt[N*N*N];
float L=1.0;
float D=0.01;
float Ux=1.0;
float Uy=0.0;
float Uz=-0.4;

float dx = L/N;
float dt = 20e-5;

float V=dx*dx*dx;
float A=dx*dx;


int main(int argc, char** argv)
{
	int iw,ie,in,is,it,ib;
	float sum=0.0;
	int row,layer;
	float diffusiondt=0.0;;
	float convectiondt=0.0;

	for (int j=0;j<N*N*N;j++)
	{
		c[j]=0.0;
		dcdt[j]=0.0;
	}
	c[N*N*N/2+N*N/2+N/2]=1.0;

	for (int i=0;i<number_of_timesteps;i++)
	{
		for (int j=0;j<N*N*N;j++)
		{
			layer=N*N*(j/(N*N));			
			row = (N*(j/N))%(N*N);
			iw=((j+(N-1))%N)+row+layer;
			ie=((j+1)%N)+row+layer;
			is=(j+N*(N-1))%(N*N)+layer;
			in=(j+N)%(N*N)+layer;	
			ib=(j+N*N*(N-1))%(N*N*N);
			it=(j+N*N)%(N*N*N);

			//cout << j << " " << iw << " "<< ie << " "<< is << " "<< in << " "<< ib << " "<< it << "   " << row << "   " << layer << endl;
	

		        diffusiondt=D*A/V/dx*(c[iw]+c[ie]+c[in]+c[is]+c[it]+c[ib]-6*c[j]);
			convectiondt=Ux*A/V*(0.5*(c[iw]+c[j])-0.5*(c[ie]+c[j]))  +Uy*A/V*(0.5*(c[is]+c[j])-0.5*(c[in]+c[j])) +Uz*A/V*(0.5*(c[ib]+c[j])-0.5*(c[it]+c[j]) );

			dcdt[j]=diffusiondt+convectiondt;	
		}

		for (int j=0;j<N*N*N;j++)
		{
			c[j]=c[j]+dcdt[j]*dt;
			sum = sum + c[j];
		}
		//cout << c[N*N*N/2+N*N/2+N/2] << " " << sum << endl;
		sum = 0.0;

	}

 	ofstream myfile;
  	myfile.open ("example.vtk");
  	myfile << "# vtk DataFile Version 3.0\n";
	myfile << "Test output\n";
	myfile << "ASCII\n";
	myfile << "DATASET RECTILINEAR_GRID\n";
	myfile << "DIMENSIONS " << N+1 << " " << N+1 << " " << N+1 << "\n";
	myfile << "X_COORDINATES " << N+1 << " float\n";
	for (int j=0;j<N+1;j++)
	{
		myfile << (j%(N+1))*L/N << "\n";
	}			
	myfile << "Y_COORDINATES " << N+1 << " float\n";
	for (int j=0;j<(N+1);j++)
	{
		myfile << (j%(N+1))*L/N << "\n";
	}
	myfile << "Z_COORDINATES " << N+1 << " float\n";
	for (int j=0;j<(N+1);j++)
	{
		myfile << (j%(N+1))*L/N << "\n";
	}

	myfile << "CELL_DATA " << N*N*N <<"\n";
	myfile << "SCALARS cellval float\n";
	myfile << "LOOKUP_TABLE default\n";
	for (int j=0;j<N*N*N;j++)
	{
		myfile << c[j] << "\n";
	}		
  	myfile.close(); 


}
