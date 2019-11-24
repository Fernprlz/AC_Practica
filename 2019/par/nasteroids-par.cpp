#include <iostream>
#include <random>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <cstdlib>

using namespace std;

const int X = 0;
const int Y = 1;

struct asteroide {
	double pX;					// Posicion eje x
	double pY;					// Posicion eje y
	double masa;
	double sig_pX;			// Siguiente pX calculada (a actualizar en pX)
	double sig_pY;			// Siguiente pY calculada (a actualizar en py)
	double vX;
	double vY;
	double sig_vX;
	double sig_vY;
};

struct planeta {
	double pX;
	double pY;
	double masa;
};

// ------------------------------------ Funciones auxiliares ------------------------------------------ //
double calcularDistanciaAsteroide(asteroide cuerpo1, asteroide cuerpo2);
double calcularPendienteAsteroide(asteroide cuerpo1, asteroide cuerpo2);
double calcularFuerzaAsteroide(asteroide cuerpo1, asteroide cuerpo2, double gravity, double distancia);
double calcularDistanciaPlaneta(planeta cuerpo1, asteroide cuerpo2);
double calcularPendientePlaneta(planeta cuerpo1, asteroide cuerpo2);
double calcularFuerzaPlaneta(planeta cuerpo1, asteroide cuerpo2, double gravity, double distancia);
void descomponerFuerzasX(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo);
void descomponerFuerzasY(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo);
void calcularNuevaAceleracion(asteroide ast, double sum_fX, double sum_fY, double *aceleracion);
void calcularNuevasVelocidades(asteroide &ast, double aceleracion[2], double time_interval);
void calcularNuevaPosicion(asteroide &ast, double time_interval);
void comprobarBordes(asteroide &ast, double width, double height);
void actualizarAsteroide(asteroide *ast, double *sum_fX, double *sum_fY);
// --------------------------------------------------------------------------------------------------- //

int main(int argc, char *argv[]) {

	//////////// V A L I D A C I O N - D E - A R G U M E N T O S ////////////
	if (argc != 5 || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 || atoi(argv[3]) < 0 || atoi(argv[4]) <= 0) {
		cout << "nasteroids-par: Wrong arguments." << endl;
		cout << "Correct use:" << endl;
		cout << "nasteroids-par num_asteroides num_iteraciones num_planetas semilla" << endl;
		return -1;
	}

	// Guardamos los argumentos en variables locales para acceder a ellas más facilmente
	int num_asteroides = atoi(argv[1]);
	int num_iteraciones = atoi(argv[2]);
	int num_planetas = atoi(argv[3]);
	int semilla = atoi(argv[4]);

	// Constantes proporcionadas por el enunciado
	const double gravity = 6.674e-5;
	const double time_interval = 0.1;
	const double dmin = 5.0;
	const double width = 200;
	const double height = 200;
	const double mean = 1000;
	const double sdm = 50;

	// Generación aleatoria del tablero utilizando la semilla
	default_random_engine random{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width, std::numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height, std::numeric_limits<double>::max())};
	normal_distribution<double> mdist{mean, sdm};

	// Crear .txt con la configuracion inicial
	ofstream fs("init_conf.txt");

	// Pasar al fichero los argumentos de la configuracion inicial
	for (int ii = 0; ii < argc - 1; ii++) {

		if (ii != 3) fs << argv[ii + 1] << " ";
		else fs << argv[ii + 1] << endl;
	}

	// Fijar precisión de 3 decimales en los archivos
	fs << setprecision(3) << fixed;

	// Vector que contendra todos los asteroides
	vector<asteroide> asteroides(num_asteroides);
	// Vector que contendra todos los planetas
	vector<planeta> planetas(num_planetas);

	// Generar asteroides llenando su vector y escribir sus parametros en el .txt con la configuracion inicial
	for(unsigned int ast = 0; ast < asteroides.size(); ast++) {
		asteroides[ast] = {xdist(random), ydist(random), mdist(random), 0, 0, 0, 0, 0, 0};
		fs << asteroides[ast].pX << " " << asteroides[ast].pY << " " << asteroides[ast].masa << endl;
	}


	// Generar planetas llenando su vector y escribir sus parametros en el .txt de la configuracion inicial
	for(unsigned int pla = 0; pla < planetas.size(); pla++) {

		if (pla % 4 == 0) {
			planetas[pla] = {0.0, ydist(random), mdist(random) * 10};
		}
		else if (pla % 4 == 1) {
			planetas[pla] = {xdist(random), 0.0, mdist(random) * 10};
		}
		else if (pla % 4 == 2) {
			planetas[pla] = {width, ydist(random), mdist(random) * 10};
		}
		else if (pla % 4 == 3) {
			planetas[pla] = {xdist(random), height, mdist(random) * 10};
		}
		fs << planetas[pla].pX << " " << planetas[pla].pY << " " << planetas[pla].masa << endl;
	}

	// Cerrar init_conf.txt
	fs.close();

	// Variables para los calculos
	double distancia = 0;
	double pendiente = 0;
	double angulo = 0;
	double aceleracion[2] = {0}; 	// eje X, eje Y
	double aux_vX = 0;						// Utilizada para intercambiar velocidades
	double aux_vY = 0;						// Utilizada para intercambiar velocidades
	double fuerza = 0;

	// Vectores que almacenan los sumatorios de fuerzas de todos los asteroides (en sus correspondientes indices)
	vector<double> sum_fX(num_asteroides);
	vector<double> sum_fY(num_asteroides);

	// Arrays de punteros de dos dimensiones. Guardan en el indice [m][n] la fuerza obtenida de la interaccion
	// entre el cuerpo m y el cuerpo n (en el correspondiente eje).
	// Los arrays serán compartidos por todos los threads. No se dan condiciones
	// de carrera porque ningun hilo tendra una pareja de indices identica
	double **sum_total_fX = new double*[num_asteroides];
	double **sum_total_fY = new double*[num_asteroides];

	// Array bidimensional utilizado para recordar las choques entre asteroides. Un 1 en choques[m][n]
	// significa que el asteroide m choca con el asteroide n.
	// Se usa para intercambiar las velocidades en el orden correcto, despues de todos los calculos
	double **choques = new double*[num_asteroides];

	// Segundas dimensiones de los arrays anteriores
	for(int ast = 0; ast < num_asteroides; ast++) {
		sum_total_fX[ast] = new double[num_asteroides];
		sum_total_fY[ast] = new double[num_asteroides];
		choques[ast] = new double[num_asteroides];
	}

	// Inicializar paralelamente los arrays anteriores
	#pragma omp parallel for
	for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {
		for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {
			sum_total_fX[ast_1][ast_2] = 0;
			sum_total_fY[ast_1][ast_2] = 0;
			choques[ast_1][ast_2] = 0;
		}
	}

	unsigned int ast_1 = 0, ast_2 = 0, pla = 0;


	// ----- Interacción asteroides-asteroides ----- //
	// Primer for loop: iteraciones introducida por el usuario
	for (int iteracion = 0; iteracion < num_iteraciones; iteracion++) {
		// Cada thread lo ejecutara tantas veces como el usuario haya introducido.
		// Las variables para los calculos se hacen privadas para cada thread, para evitar condiciones de carrera.
		#pragma omp parallel for private(ast_2, fuerza, distancia, pendiente, angulo)
		// Segundo for loop: asteroides
		for (ast_1 = 0; ast_1 < asteroides.size(); ast_1++) {

			// Tercer for loop: interaccion entre el asteroide ast_1 y todos los siguientes
			// Con esta se comprueba la interaccion entre todos los asteroides
			for (ast_2 = ast_1 + 1; ast_2 < asteroides.size(); ast_2++) {

				// Distancia entre asteroides
				distancia = calcularDistanciaAsteroide(asteroides[ast_1], asteroides[ast_2]);

				// Si distancia <= dmin (5), se señala el choque entre ast_1 y ast_2 para posteriormente intercambiar sus velocidades
				if (distancia <= dmin) {
					if (ast_1 > 0) {
						choques[ast_1][ast_2] = 1;
					}

				// Caso de que la distancia entre asteroides sea mayor que 5:
				// calcular pendiente, angulo y fuerza
				} else {
					// Calcular pendiente. Si es mayor que 1, se trunca a 1. Si es menor que -1, se trunca a -1
					pendiente =  calcularPendienteAsteroide(asteroides[ast_1], asteroides[ast_2]);

					// Calcular angulo
					angulo = atan(pendiente);

					// Calcular fuerza. Si es mayor que 100, se trunca a este valor
					fuerza = calcularFuerzaAsteroide(asteroides[ast_1], asteroides[ast_2], gravity, distancia);

					// Descomponer fuerzas: eje X / eje Y
					descomponerFuerzasX(&sum_total_fX[ast_1][ast_2], &sum_total_fX[ast_2][ast_1], fuerza, angulo);
					descomponerFuerzasY(&sum_total_fY[ast_1][ast_2], &sum_total_fY[ast_2][ast_1], fuerza, angulo);
				}
			}
		}

		// Copiar de forma paralela a sum_fX[ast_1] y sum_fY[ast_1] el sumatorio de todas las fuerzas
		// del asteroide correspondiente
		#pragma omp parallel for
		for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {
			for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {
				sum_fX[ast_1] += sum_total_fX[ast_1][ast_2];
				sum_fY[ast_1] += sum_total_fY[ast_1][ast_2];
			}
		}

		// Nested loop que intercambia las velocidades de los asteroides que hayan chocado.
		// Para ello se utiliza la matriz de choques: si se encuentra un 1 en choques[m][n],
		// se intercambian las velocidades de los asteroides m y n.
		// Se ejecuta de forma secuencial porque accede a las velocidades de los asteroides.
		for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {
			for (int ast_2 = ast_1 + 1; ast_2 < num_asteroides; ast_2++) {
				if (choques[ast_1][ast_2] == 1) {
					aux_vX = asteroides[ast_1].vX;
					aux_vY = asteroides[ast_1].vY;
					asteroides[ast_1].vX = asteroides[ast_2].vX;
					asteroides[ast_1].vY = asteroides[ast_2].vY;
					asteroides[ast_2].vX = aux_vX;
					asteroides[ast_2].vY = aux_vY;
				}
			}
		}


		// ----- Interacción planetas-asteroides ----- //
		#pragma omp parallel for private(pla, fuerza, distancia, pendiente, angulo)
		// Evaluar interaccion entre cada asteroide con cada planeta
		for (unsigned int ast = 0; ast < asteroides.size(); ast++) {
			for (pla = 0; pla < planetas.size(); pla++) {

				// Calcular distancia entre planeta y asteroide
				distancia = calcularDistanciaPlaneta(planetas[pla], asteroides[ast]);
				// Calcular pendiente. Si es mayor que 1, se trunca a 1. Si es menor que -1, se trunca a -1
			  pendiente = calcularPendientePlaneta(planetas[pla], asteroides[ast]);

				// Calcular angulo
				angulo = atan(pendiente);

				// Calcular fuerza. Si es mayor que 100, se trunca a este valor
				fuerza = calcularFuerzaPlaneta(planetas[pla], asteroides[ast], gravity, distancia);

				// Sumar cada componente de la fuerza al sumatorio de cada eje (X, Y)
				sum_fX[ast] += fuerza * cos(angulo);
				sum_fY[ast] += fuerza * sin(angulo);
			}
		}

		// Calcular nuevos parametros de los asteroides de forma paralela
		// Cada thread tiene su copia del array aceleracion, que contiene el eje X y el Y
		#pragma omp parallel for private(aceleracion)
		for (unsigned int ast = 0; ast < asteroides.size(); ast++) {

			// Nueva aceleracion
			calcularNuevaAceleracion(asteroides[ast], sum_fX[ast], sum_fY[ast], aceleracion);

			// Nueva velocidad
			calcularNuevasVelocidades(asteroides[ast], aceleracion, time_interval);

			// Nueva posicion
			calcularNuevaPosicion(asteroides[ast], time_interval);

			// Si se da rebote contra un borde, se coloca el asteroide a 5 unidades de distancia y se cambia su velocidad
			comprobarBordes(asteroides[ast], width, height);
		}

		// Actualizar paralelamente los asteroides con las nuevas velocidades. Fuerzas reseteadas a 0
		#pragma omp parallel for
		for (unsigned int ast_1 = 0; ast_1 < asteroides.size(); ast_1++) {
			actualizarAsteroide(&asteroides[ast_1], &sum_fX[ast_1], &sum_fY[ast_1]);
			for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {
				sum_total_fX[ast_1][ast_2] = 0;
				sum_total_fY[ast_1][ast_2] = 0;
				choques[ast_1][ast_2] = 0;
			}
		}
	}

	// Crear fichero para output (out.txt)
	// Fijar precision en 3 decimales
	ofstream fs2("out.txt");
	fs2 << setprecision(3) << fixed;

	// Escribir en out.txt los parametros finales de todos los asteroides
	for (unsigned int i = 0; i < asteroides.size(); i++) {

		fs2 << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].vX << " " << asteroides[i].vY
		<< " " << asteroides[i].masa << endl;
	}

	// Cerrar out.txt
	fs2.close();

	// Liberar memoria
	#pragma omp parallel for
	for(int i = 0 ; i < num_asteroides; i++) {

		delete[] sum_total_fX[i];
		delete[] sum_total_fY[i];
		delete[] choques[i];
	}

	delete[] sum_total_fX;
	delete[] sum_total_fY;
	delete[] choques;

	return 0;
}




/////////////////////////////// Cuerpo de las funciones auxiliares ///////////////////////////////////
double calcularDistanciaAsteroide(asteroide cuerpo1, asteroide cuerpo2){
	return sqrt(pow(cuerpo1.pX - cuerpo2.pX, 2.0) + pow(cuerpo1.pY - cuerpo2.pY, 2.0));
}

double calcularPendienteAsteroide(asteroide cuerpo1, asteroide cuerpo2){
	double pendiente = (cuerpo1.pY - cuerpo2.pY) / (cuerpo1.pX - cuerpo2.pX);

	// Si pendiente > 1, se trunca a 1
	if (pendiente > 1){
		pendiente = 1;
	}
	// Si pendiente < -1, se trunca a -1
	else if (pendiente < -1) {
		pendiente = -1;
	}

	return pendiente;
}

double calcularFuerzaAsteroide(asteroide cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
	// Si fuerza > 100, se trunca a 100
	if (fuerza > 100){
		fuerza = 100;
	}
	return fuerza;
}

double calcularDistanciaPlaneta(planeta cuerpo1, asteroide cuerpo2){
	return sqrt(pow(cuerpo1.pX - cuerpo2.pX, 2.0) + pow(cuerpo1.pY - cuerpo2.pY, 2.0));

}

double calcularPendientePlaneta(planeta cuerpo1, asteroide cuerpo2){
	double pendiente = (cuerpo1.pY - cuerpo2.pY) / (cuerpo1.pX - cuerpo2.pX);

	// Si pendiente > 1, se trunca a 1
	if (pendiente > 1){
		pendiente = 1;
	}
	// Si pendiente < -1, se trunca a -1
	else if (pendiente < -1) {
		pendiente = -1;
	}

	return pendiente;
}
double calcularFuerzaPlaneta(planeta cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
	// Si fuerza > 100, se trunca a 100
		if (fuerza > 100){
			fuerza = 100;
		}
		return fuerza;
}

void descomponerFuerzasX(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo){
	*sum_total_f1 += fuerza * cos(angulo);
	*sum_total_f2 -= fuerza * cos(angulo);
}

void descomponerFuerzasY(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo){
	*sum_total_f1 += fuerza * sin(angulo);
	*sum_total_f2 -= fuerza * sin(angulo);
}

void calcularNuevaAceleracion(asteroide ast, double sum_fX, double sum_fY, double *aceleracion){
	aceleracion[X] = (1/ast.masa) * sum_fX;
	aceleracion[Y] = (1/ast.masa) * sum_fY;
}
void calcularNuevasVelocidades(asteroide &ast, double aceleracion[2], double time_interval){
	ast.sig_vX = ast.vX + (aceleracion[X] * time_interval);
	ast.sig_vY = ast.vY + (aceleracion[Y] * time_interval);
}
void calcularNuevaPosicion(asteroide &ast, double time_interval){
	ast.sig_pX = ast.pX + (ast.sig_vX * time_interval);
	ast.sig_pY = ast.pY + (ast.sig_vY * time_interval);
}
void comprobarBordes(asteroide &ast, double width, double height){
	if (ast.sig_pX <= 0) {
		ast.sig_pX = 5;
		ast.sig_vX *= (-1);
	}
	if (ast.sig_pX >= width) {
		ast.sig_pX = width - 5;
		ast.sig_vX *= (-1);
	}
	if (ast.sig_pY <= 0) {
		ast.sig_pY = 5;
		ast.sig_vY *= (-1);
	}
	if (ast.sig_pY >= height) {
		ast.sig_pY = height - 5;
		ast.sig_vY *= (-1);
	}
}
void actualizarAsteroide(asteroide *ast, double *sum_fX, double *sum_fY){
	ast -> pX = ast -> sig_pX;
	ast -> pY = ast -> sig_pY;
	ast -> vX = ast -> sig_vX;
	ast -> vY = ast -> sig_vY;
	*sum_fX = 0;
	*sum_fY = 0;
}
