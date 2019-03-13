#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <algorithm>
#include <memory>
#include <math.h> 

using namespace std;

void Solver::Init(unsigned N, float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;
	this->N = N;
}

void Solver::FreeData(void)
{
	//TODO: Libera los buffers de memoria.

	if (u) free(u);
	if (v) free(v);
	if (dens) free(dens);
	if (u_prev) free(u_prev);
	if (v_prev) free(v_prev);
	if (dens_prev) free(dens_prev);
}

void Solver::ClearData(void)
{

	//TODO: Borra todo el contenido de los buffers
	memset(u, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(v, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(dens, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(u_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(v_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(dens_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
}

bool Solver::AllocateData(void)
{
	//TODO:
	//Reservamos memoria, en caso de fallo devlvemos false.
	//Antes de devolver true, hay que limpiar la memoria reservada con un ClearData().
	try
	{
		u = (float*)malloc((N + 2)*(N + 2) * sizeof(float));
		v = (float*)malloc((N + 2)*(N + 2) * sizeof(float));
		dens = (float*)malloc((N + 2)*(N + 2) * sizeof(float));
		u_prev = (float*)malloc((N + 2)*(N + 2) * sizeof(float));
		v_prev = (float*)malloc((N + 2)*(N + 2) * sizeof(float));
		dens_prev = (float*)malloc((N + 2)*(N + 2) * sizeof(float));

		ClearData();
		return true;
	}
	catch (const std::exception&)
	{
		printf("Fallo al asignar memoria");
		return false;
	}
}

void Solver::ClearPrevData()
{
	//TODO: Borra el contenido de los buffers _prev
	memset(u_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(v_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
	memset(dens_prev, 0.0f, (N + 2)*(N + 2) * sizeof(float));
}

void Solver::AddDensity(unsigned x, unsigned y, float source)
{
	//TODO: Añade el valor de source al array de densidades. Sería interesante usar la macro: XY_TO_ARRAY
	dens_prev[XY_TO_ARRAY(x, y)] = source;
}

void Solver::AddVelocity(unsigned x, unsigned y, float forceX, float forceY)
{
	//TODO: Añade el valor de fuerza a sus respectivos arrays. Sería interesante usar la macro: XY_TO_ARRAY
	u_prev[XY_TO_ARRAY(x, y)] = forceX;
	v_prev[XY_TO_ARRAY(x, y)] = forceY;
}

void Solver::Solve()
{
	VelStep();
	DensStep();
}

void Solver::DensStep()
{

	AddSource(dens, dens_prev);			//Adding input density (dens_prev) to final density (dens).
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Diffuse(0, dens, dens_prev);		//Writing result in dens because we made the swap before. bi = dens_prev. The initial trash in dens matrix, doesnt matter, because it converges anyways.
	SWAP(dens_prev, dens)				//Swapping matrixes, because we want save the next result in dens, not in dens_prev.
	Advect(0, dens, dens_prev, u, v);	//Advect phase, result in dens.
}

void Solver::VelStep()
{
	if (soplador) {
		u[XY_TO_ARRAY(1, N / 2)] = 1;
	}

	AddSource(u, u_prev);
	AddSource(v, v_prev);
	SWAP(u_prev, u)
	SWAP(v_prev, v)
	Diffuse(1, u, u_prev);
	Diffuse(2, v, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
	SWAP(u_prev, u)
	SWAP(v_prev, v)
	Advect(1, u, u_prev, u_prev, v_prev);
	Advect(2, v, v_prev, u_prev, v_prev);
	Project(u, v, u_prev, v_prev);		//Mass conserving.
}

void Solver::AddSource(float * base, float * source)
{
	//TODO: Teniendo en cuenta dt (Delta Time), incrementar el array base con nuestro source. 
	//      Esto sirve tanto para añadir las nuevas densidades como las nuevas fuerzas.
	FOR_EACH_CELL
		base[XY_TO_ARRAY(i, j)] += source[XY_TO_ARRAY(i, j)] * dt;
	END_FOR
}


void Solver::SetBounds(int b, float * x)
{
	/*TODO:
	Input b: 0, 1 or 2.
		0: borders = same value than the inner value.
		1: x axis borders inverted, y axis equal.
		2: y axis borders inverted, x axis equal.
		Corner values allways are mean value between associated edges.
	*/
	int i, j;
	FOR_EACH_CELL
			if (i == 0 && j < N + 1)
			{
				switch (b) {
				case 0:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i + 1, j)];
					break;
				case 1:
					x[XY_TO_ARRAY(i, j)] = -x[XY_TO_ARRAY(i + 1, j)];
					break;
				case 2:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i + 1, j)];
					break;
				}

			}
			else if (j == 0 && i < N + 1)
			{
				switch (b) {
				case 0:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i, j + 1)];
					break;
				case 1:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i, j + 1)];
					break;
				case 2:
					x[XY_TO_ARRAY(i, j)] = -x[XY_TO_ARRAY(i, j + 1)];
					break;
				}
			}
			else if (i == N + 1 && j < N + 1)
			{
				switch (b) {
				case 0:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i - 1, j)];
					break;
				case 1:
					x[XY_TO_ARRAY(i, j)] = -x[XY_TO_ARRAY(i - 1, j)];
					break;
				case 2:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i - 1, j)];
					break;
				}
			}
			else if (j == N + 1 && i < N + 1)
			{
				switch (b) {
				case 0:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i, j - 1)];
					break;
				case 1:
					x[XY_TO_ARRAY(i, j)] = x[XY_TO_ARRAY(i, j - 1)];
					break;
				case 2:
					x[XY_TO_ARRAY(i, j)] = -x[XY_TO_ARRAY(i, j - 1)];
					break;
				}
			}
	END_FOR

	x[XY_TO_ARRAY(0, 0)] = (x[XY_TO_ARRAY(1, 0)] + x[XY_TO_ARRAY(0, 1)]) / 2.0;

	x[XY_TO_ARRAY(0, N + 1)] = (x[XY_TO_ARRAY(0, N)] + x[XY_TO_ARRAY(1, N + 1)]) / 2.0;

	x[XY_TO_ARRAY(N + 1, 0)] = (x[XY_TO_ARRAY(N, 0)] + x[XY_TO_ARRAY(N + 1, 1)]) / 2.0;

	x[XY_TO_ARRAY(N + 1, N + 1)] = (x[XY_TO_ARRAY(N, N + 1)] + x[XY_TO_ARRAY(N + 1, N)]) / 2.0;

}

/*
Gauss Seidel -> Matrix x and x0
*/
void Solver::LinSolve(int b, float * x, float * x0, float aij, float aii) {
	//TODO: Se recomienda usar FOR_EACH_CELL, END_FOR y XY_TO_ARRAY.
	for (int gs = 0; gs < 20; ++gs) {
		FOR_EACH_CELL
			if (jacobi) {
				x[XY_TO_ARRAY(i, j)] = 
				(-aij * (-x0[XY_TO_ARRAY(i, j - 1)] 
						- x0[XY_TO_ARRAY(i - 1, j)] 
						- x0[XY_TO_ARRAY(i + 1, j)] 
						- x0[XY_TO_ARRAY(i, j + 1)]) 
						+ x0[XY_TO_ARRAY(i, j)]) / aii;
			}
			else{
				x[XY_TO_ARRAY(i, j)] =
				(-aij * (-x[XY_TO_ARRAY(i, j - 1)] //parte del sumatorio de la funcion de gauss seidel
						- x[XY_TO_ARRAY(i - 1, j)]
						- x[XY_TO_ARRAY(i + 1, j)]
						- x[XY_TO_ARRAY(i, j + 1)])
						+ x0[XY_TO_ARRAY(i, j)]) / aii; // sumamos el valor previo y lo dividimos entre aii .
			}
			
		END_FOR
		SetBounds(b, x); //hacemos el setbound en cada iteracion de las 20. consejo de jesus
		memccpy(x0, x, (N + 2) * (N + 2), sizeof(float)); //cambiamos el buffer para dejar el x0 como el prev_buffer para la siguiente iteracion
	}
}


/*
Nuestra función de difusión solo debe resolver el sistema de ecuaciones simplificado a las celdas contiguas de la casilla que queremos resolver,
por lo que solo con la entrada de dos valores, debemos poder obtener el resultado.
*/
void Solver::Diffuse(int b, float * x, float * x0) {
	//TODO: Solo necesitaremos pasar dos parámetros a nuestro resolutor de sistemas de ecuaciones de Gauss Seidel. Calculamos dichos valores y llamamos a la resolución del sistema.
	float aij = diff * dt * N * N;
	float aii = 1 + 4 * aij;
	LinSolve(b, x, x0, aij, aii);
}

/*
d is overwrited with the initial d0 data and affected by the u & v vectorfield.
Hay que tener en cuenta que el centro de las casillas representa la posición entera dentro de la casilla, por lo que los bordes estan
en las posiciones x,5.
*/
void Solver::Advect(int b, float * d, float * d0, float * u, float * v) {
	//TODO: Se aplica el campo vectorial realizando una interploación lineal entre las 4 casillas más cercanas donde caiga el nuevo valor.

	float director, deltaX, deltaY, oX_f, oY_f;

	director = -dt * N;
		int oX, oY;

	FOR_EACH_CELL
	oX_f = i + (u[XY_TO_ARRAY(i, j)] * director);
	oY_f = j + (v[XY_TO_ARRAY(i, j)] * director);
	oX = (int)oX_f;
	oY = (int)oY_f;

	if (oX > 0 && oX < N+1 && oY > 0 && oY < N+1) {

		deltaX = fabs(oX_f - oX);
		deltaY = fabs(oY_f - oY);
		d[XY_TO_ARRAY(i, j)] = d0[XY_TO_ARRAY(oX, oY)] * (1 - deltaX) * (1 - deltaY) +
			d0[XY_TO_ARRAY(oX + 1, oY)] * deltaX * (1 - deltaY) +
			d0[XY_TO_ARRAY(oX, oY + 1)] * (1 - deltaX) * deltaY +
			d0[XY_TO_ARRAY(oX + 1, oY + 1)] * deltaX * deltaY;
	}
	else {
		d[XY_TO_ARRAY(i, j)] = 0;
	}
	END_FOR
		SetBounds(b, d);
}

/*
Se encarga de estabilizar el fluido y hacerlo conservativo de masa. Se usa solo en las matrices de velocidades.
No necesaria implementación por su complejidad.
*/
void Solver::Project(float * u, float * v, float * p, float * div)
{

	FOR_EACH_CELL
		div[XY_TO_ARRAY(i, j)] = -0.5f*(u[XY_TO_ARRAY(i + 1, j)] - u[XY_TO_ARRAY(i - 1, j)] + v[XY_TO_ARRAY(i, j + 1)] - v[XY_TO_ARRAY(i, j - 1)]) / N;
	p[XY_TO_ARRAY(i, j)] = 0;
	END_FOR
		SetBounds(0, div);
	SetBounds(0, p);

	LinSolve(0, p, div, 1, 4);

	//Aproximamos: Laplaciano de q a su gradiente.
	FOR_EACH_CELL
		u[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i + 1, j)] - p[XY_TO_ARRAY(i - 1, j)]);
	v[XY_TO_ARRAY(i, j)] -= 0.5f*N*(p[XY_TO_ARRAY(i, j + 1)] - p[XY_TO_ARRAY(i, j - 1)]);
	END_FOR
		SetBounds(1, u);
	SetBounds(2, v);
}