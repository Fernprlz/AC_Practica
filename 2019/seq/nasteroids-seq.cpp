#include <iostream>
#include <random>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

const int X = 0;
const int Y = 1;

struct asteroide {
	double pX;								// Posicion eje x
	double pY;								// Posicion eje y
	double masa;
	double sig_pX;						// Siguiente pX calculada (a actualizar en pX)
	double sig_pY;						// Siguiente pY calculada (a actualizar en py)
	double sum_fX;
	double sum_fY;
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
void actualizarAsteroide(asteroide *ast);
double calcularDistanciaAsteroide(asteroide cuerpo1, asteroide cuerpo2);
double calcularDistanciaPlaneta(planeta cuerpo1, asteroide cuerpo2);
double calcularPendienteAsteroide(asteroide cuerpo1, asteroide cuerpo2);
double calcularPendientePlaneta(planeta cuerpo1, asteroide cuerpo2);
double calcularFuerzaAsteroide(asteroide cuerpo1, asteroide cuerpo2, double gravity, double distancia);
double calcularFuerzaPlaneta(planeta cuerpo1, asteroide cuerpo2, double gravity, double distancia);
void descomponerFuerzas(asteroide &ast1, asteroide &ast2, double fuerza, double angulo);
void calcularNuevaAceleracion(asteroide ast, double *aceleracion);
void calcularNuevasVelocidades(asteroide &ast, double aceleracion[2], double time_interval);
void calcularNuevaPosicion(asteroide &ast, double time_interval);
void comprobarBordes(asteroide &ast, double width, double height);
// --------------------------------------------------------------------------------------------------- //

int main(int argc, char *argv[]) {

	//////////// V A L I D A C I O N - D E - A R G U M E N T O S ////////////
	if (argc != 5 || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 || atoi(argv[3]) < 0 || atoi(argv[4]) <= 0) {
		cout << "nasteroids-seq: Wrong arguments." << endl;
		cout << "Correct use:" << endl;
		cout << "nasteroids-seq num_asteroides num_iteraciones num_planetas semilla" << endl;
		return -1;
	}

	// Guardamos los argumentos en variables locales para acceder a ellas más facilmente
	int num_asteroides = atoi(argv[1]);
	int num_iteraciones = atoi(argv[2]);
	int num_planetas = atoi(argv[3]);
	int semilla = atoi(argv[4]);

	// Constantes proporcionadas por el enunciado -
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
 	for (int i = 0; i < argc - 1; i++) {
 		if (i != 3){
			fs << argv[i+1] << " ";
		} else {
			fs << argv[i+1] << endl;
		}
 	}

 	// Fijar precisión de 3 decimales en los archivos
 	fs << setprecision(3) << fixed;

	// Vector que contendra todos los asteroides
	vector<asteroide> asteroides(num_asteroides);
	// Vector que contendra todos los planetas
	vector<planeta> planetas(num_planetas);

	// Generar asteroides llenando su vector y escribir sus parametros en el .txt con la configuracion inicial
	for(unsigned int i = 0; i < asteroides.size(); i++) {
		asteroides[i] = {xdist(random), ydist(random), mdist(random), 0, 0, 0, 0, 0, 0, 0, 0};
		fs << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].masa << endl;
	}

	// Generar planetas llenando su vector y escribir sus parametros en el .txt de la configuracion inicial
	for(unsigned int i = 0; i < planetas.size(); i++) {

		if (i % 4 == 0) {
			planetas[i] = {0.0, ydist(random), mdist(random) * 10};
		}
		else if (i % 4 == 1) {
			planetas[i] = {xdist(random), 0.0, mdist(random) * 10};
		}
		else if (i % 4 == 2) {
			planetas[i] = {200.0, ydist(random), mdist(random) * 10};
		}
		else if (i % 4 == 3) {
			planetas[i] = {xdist(random), 200.0, mdist(random) * 10};
		}
		fs << planetas[i].pX << " " << planetas[i].pY << " " << planetas[i].masa << endl;
	}

	// Cerrar init_conf.txt */
	fs.close();

	// Variables para los calculos
	double distancia = 0;
	double pendiente = 0;
	double angulo = 0;
	double aceleracion[2] = {0}; 	// eje X, eje Y
	double vX_aux = 0;						// Utilizada para intercambiar velocidades
	double vY_aux = 0;						// Utilizada para intercambiar velocidades
	double fuerza = 0;

	// Primer for loop: iteraciones introducida por el usuario
	for (int i = 0; i < num_iteraciones; i++) {

		// Segundo for loop: asteroides
		for (unsigned int aa = 0; aa < asteroides.size(); aa++) {

			// Tercer for loop: interaccion entre el asteroide aa y todos los siguientes
			// Con esta se comprueba la interaccion entre todos los asteroides
			for (unsigned int bb = aa+1; bb < asteroides.size(); bb++) {

				// Distancia entre asteroides
				distancia = calcularDistanciaAsteroide(asteroides[aa], asteroides[bb]);

				// Caso de que la distancia entre asteroides sea menor que 5:
				// se intercambian las velocidades de los dos asteroides
				if (distancia <= dmin && aa > 0) {
						vX_aux = asteroides[aa].vX;
						vY_aux = asteroides[aa].vY;
						asteroides[aa].vX = asteroides[bb].vX;
						asteroides[aa].vY = asteroides[bb].vY;
						asteroides[bb].vX = vX_aux;
						asteroides[bb].vY = vY_aux;
				}

				// Caso de que la distancia entre asteroides sea mayor que 5:
				// calcular pendiente, angulo y fuerza
				else {

					// Calcular pendiente. Si es mayor que 1, se trunca a 1. Si es menor que -1, se trunca a -1
					pendiente = calcularPendienteAsteroide(asteroides[aa], asteroides[bb]);

					// Calcular angulo
					angulo = atan(pendiente);

					// Calcular fuerza. Si es mayor que 100, se trunca a este valor
					fuerza = calcularFuerzaAsteroide(asteroides[aa], asteroides[bb], gravity, distancia);

					// Descomponer fuerzas: eje X / eje Y
					descomponerFuerzas(asteroides[aa], asteroides[bb], fuerza, angulo);
				}
			}
		}

		// Interaccion de planetas (pp) con asteroides (aa)
		for (unsigned int pp = 0; pp < planetas.size(); pp++) {
			for (unsigned int aa = 0; aa < asteroides.size(); aa++) {

				// Calcular distancia entre planeta y asteroide
				distancia = calcularDistanciaPlaneta(planetas[pp], asteroides[aa]);

				// Calcular pendiente. Si es mayor que 1, se trunca a 1. Si es menor que -1, se trunca a -1
				pendiente = calcularPendientePlaneta(planetas[pp], asteroides[aa]);

				// Calcular angulo
				angulo = atan(pendiente);

				// Calcular fuerza. Si es mayor que 100, se trunca a este valor
				fuerza = calcularFuerzaPlaneta(planetas[pp], asteroides[aa], gravity, distancia);

				// Sumar cada componente de la fuerza al sumatorio de cada eje (X, Y)
				asteroides[aa].sum_fX += fuerza * cos(angulo);
				asteroides[aa].sum_fY += fuerza * sin(angulo);
			}
		}

		// Calcular nuevos parametros de los asteroides
		for (unsigned int aa = 0; aa < asteroides.size(); aa++) {

			// Nueva aceleracion
			calcularNuevaAceleracion(asteroides[aa], aceleracion);

			// Nueva velocidad
			calcularNuevasVelocidades(asteroides[aa], aceleracion, time_interval);

			// Nueva posicion
			calcularNuevaPosicion(asteroides[aa], time_interval);

			// Si se da rebote contra un borde, se coloca el asteroide a 5 unidades de distancia y se cambia su velocidad
			comprobarBordes(asteroides[aa], width, height);
		}

		// Actualizar asteroides con las nuevas velocidades. Fuerzas reseteadas a 0
		for (unsigned int aa = 0; aa< asteroides.size(); aa++) {
			actualizarAsteroide(&asteroides[aa]);
			cout << "Velocidad del asteroide [" << aa << "]: " << asteroides[aa].sig_vX << ", " << asteroides[aa].sig_vY << "\n";
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

	return 0;
}





/////////////////////////////// Cuerpo de las funciones auxiliares ///////////////////////////////////
void actualizarAsteroide(asteroide *ast){
	ast -> pX = ast -> sig_pX;
	ast -> pY = ast -> sig_pY;
	ast -> vX = ast -> sig_vX;
	ast -> vY = ast -> sig_vY;
	ast -> sum_fX = 0;
	ast -> sum_fY = 0;
}

double calcularDistanciaAsteroide(asteroide cuerpo1, asteroide cuerpo2){
	return sqrt(pow(cuerpo1.pX - cuerpo2.pX, 2.0) + pow(cuerpo1.pY - cuerpo2.pY, 2.0));
}

double calcularDistanciaPlaneta(planeta cuerpo1, asteroide cuerpo2){
	return sqrt(pow(cuerpo1.pX - cuerpo2.pX, 2.0) + pow(cuerpo1.pY - cuerpo2.pY, 2.0));
}

double calcularPendienteAsteroide(asteroide cuerpo1, asteroide cuerpo2){

	double pendiente = (cuerpo1.pY - cuerpo2.pY) / (cuerpo1.pX - cuerpo2.pX);

	/* Si la pendiente es mayor a 1, se fija su valor a 1*/
	if (pendiente > 1){
		pendiente = 1;
	}
	/* Si la pendiente es menor a -1, se fija su valor a -1*/
	else if (pendiente < -1) {
		pendiente = -1;
	}

	return pendiente;
}

double calcularPendientePlaneta(planeta cuerpo1, asteroide cuerpo2){

	double pendiente = (cuerpo1.pY - cuerpo2.pY) / (cuerpo1.pX - cuerpo2.pX);

	/* Si la pendiente es mayor a 1, se fija su valor a 1*/
	if (pendiente > 1){
		pendiente = 1;
	}
	/* Si la pendiente es menor a -1, se fija su valor a -1*/
	else if (pendiente < -1) {
		pendiente = -1;
	}

	return pendiente;
}

double calcularFuerzaAsteroide(asteroide cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
	/* Si la fuerza resultante es mayor a 100, se trunca a este valor */
	if (fuerza > 100){
		fuerza = 100;
	}
	return fuerza;
}

double calcularFuerzaPlaneta(planeta cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
	/* Si la fuerza resultante es mayor a 100, se trunca a este valor */
	if (fuerza > 100){
		fuerza = 100;
	}
	return fuerza;
}

void descomponerFuerzas(asteroide &ast1, asteroide &ast2, double fuerza, double angulo){
		/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje x */
	ast1.sum_fX += fuerza * cos(angulo);
	ast2.sum_fX -= fuerza * cos(angulo);

	/* La fuerza resultante por el seno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje y */
	ast1.sum_fY += fuerza * sin(angulo);
	ast2.sum_fY -= fuerza * sin(angulo);
}

void calcularNuevaAceleracion(asteroide ast, double *aceleracion){
	aceleracion[X] = (1/ast.masa) * ast.sum_fX;
	aceleracion[Y] = (1/ast.masa) * ast.sum_fY;
}

// Paso el array de aceleracion con el tamaño especificado por rollos de backward compatibility con c
// que hacen que se pierda la informacion del tamaño de la array:
// https://stackoverflow.com/questions/6165449/c-error-invalid-types-intint-for-array-subscript
void calcularNuevasVelocidades(asteroide &ast, double aceleracion[2], double time_interval) {
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
