#include <iostream>
#include <random>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

// Constantes para programar bien porque molo - Linares
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

// ------ Funciones auxiliares ------ //
void actualizarAsteroide(asteroide ast);
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

int main(int argc, char *argv[]) {

	/* Comprueba el numero de argumentos e imprime el mensaje de error correspondiente */
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

	// Constantes proporcionadas por el enunciado ->> NO LAS ASEMOS GLOVALES PK SE FUCKEA EL PARALELO - GOGLE
	const double gravity = 6.674e-5;
	const double time_interval = 0.1;
	const double dmin = 5.0;
	const double width = 200;
	const double height = 200;
	const double mean = 1000;
	const double sdm = 50;

	// Generación aleatoria del tablero utilizando la semilla
	default_random_engine re{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width, std::numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height, std::numeric_limits<double>::max())};
	normal_distribution<double> mdist{mean, sdm};

	/* Crear fichero de salida */
   	ofstream fs("init_conf.txt");

   	/* Escribir argumentos en el fichero init_conf.txt */
   	for (int i = 0; i < argc - 1; i++) {

   		if (i != 3) fs << argv[i+1] << " ";
   		else fs << argv[i+1] << endl;
   	}

   	/* Fija 3 decimales al escribir en los ficheros */
   	fs << setprecision(3) << fixed;

		// Vector que almacena los asteroides
	vector<asteroide> asteroides(num_asteroides);
	vector<planeta> planetas(num_planetas); /* Vector que almacena los planetas */

	/* Generacion de los asteroides y escritura de sus datos en el fichero init_conf.txt */
	for(unsigned int i = 0; i < asteroides.size(); i++) {

		asteroides[i] = {xdist(re), ydist(re), mdist(re), 0, 0, 0, 0, 0, 0, 0, 0};
		fs << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].masa << endl;
	}

	/* Generacion de los planetas y escritura de sus datos en el fichero init_conf.txt */
	for(unsigned int i = 0; i < planetas.size(); i++) {

		if (i % 4 == 0) {
			planetas[i] = {0.0, ydist(re), mdist(re) * 10};
		}
		else if (i % 4 == 1) {
			planetas[i] = {xdist(re), 0.0, mdist(re) * 10};
		}
		else if (i % 4 == 2) {
			planetas[i] = {200.0, ydist(re), mdist(re) * 10};
		}
		else if (i % 4 == 3) {
			planetas[i] = {xdist(re), 200.0, mdist(re) * 10};
		}
		fs << planetas[i].pX << " " << planetas[i].pY << " " << planetas[i].masa << endl;
	}

	/* Cierra el fichero init_cong.txt */
	fs.close();

	/* Definir e inicializar variables auxiliares para los calculos */
	double distancia = 0;
	double pendiente = 0;
	double angulo = 0;
	double aceleracion[2]; // 2 componentes, X e Y
	double vX_aux = 0;
	double vY_aux = 0;
	double fuerza = 0;

	// Bucle exterior que ejecuta las iteraciones de la simulación
	for (int i = 0; i < num_iteraciones; i++) {

		// Primer bucle que itera sobre los asteroides
		for (unsigned int j = 0; j < asteroides.size(); j++) {

			/* Se evalua la interaccion entre el asteroide "j" y los demas que tengan índice mayor (a partir de j+1) */
			for (unsigned int k = j+1; k < asteroides.size(); k++) {

				/* Calcula la distancia entre el asteroide "j" y el asteroide "k" */
				distancia = calcularDistanciaAsteroide(asteroides[j], asteroides[k]);

				/* Si la distancia entre ambos asteroides es menor a 5 se intercambian las velocidades de la iteracion anterior */
				if (distancia <= dmin && j > 0) {
						vX_aux = asteroides[j].vX;
						vY_aux = asteroides[j].vY;
						asteroides[j].vX = asteroides[k].vX;
						asteroides[j].vY = asteroides[k].vY;
						asteroides[k].vX = vX_aux;
						asteroides[k].vY = vY_aux;
				}

				/* Si la distancia es mayor que 5 se toma en cuenta la fuerza de atraccion entre ellos */
				else {

					/* Calcula la pendiente entre el asteroide "j" y el asteroide "k", equilibrandola si es mayor a 1 o menor que -1 */
					pendiente = calcularPendienteAsteroide(asteroides[j], asteroides[k]);

					/* El angulo entre asteroides sera la ArcTan(pendiente) */
					angulo = atan(pendiente);

					/* Calculo de la fuerza resultante y si sobre pasa el valor 100, se trunca*/
					fuerza = calcularFuerzaAsteroide(asteroides[j], asteroides[k], gravity, distancia);

					// Se descomponen las fuerzas en componentes x e y.
					descomponerFuerzas(asteroides[j], asteroides[k], fuerza, angulo);
				}
			}
		}

		//TODO: PROBAR A CAMBIAR EL ORDEN DE LOS BUCLES (EXT TO INT, INT TO EXT), ES CONCEPTUALMENTE LO MISMO?
		for (unsigned int q = 0; q < planetas.size(); q++) {
			/* Para el mismo asteroide "z", despues de evaluar la interaccion con el resto de asteroides, se hace lo propio
			para con todos los planetas */
			for (unsigned int z = 0; z < asteroides.size(); z++) {

				/* Calcula la distancia entre el asteroide "z" y el planeta "q" */
				distancia = calcularDistanciaPlaneta(planetas[q], asteroides[z]);

				/* Calculo de la pendiente entre elementos y equilibro el rollo */
				pendiente = calcularPendientePlaneta(planetas[q], asteroides[z]);

				/* EL angulo entre asteroide y planeta sera el ArcTan(pendiente) */
				angulo = atan(pendiente);

				/* Calculo de la fuerza en el resultante entre ambos elementos */
				fuerza = calcularFuerzaPlaneta(planetas[q], asteroides[z], gravity, distancia);

				/* El producto de la fuerza por el coseno del angulo se agrega positvamente al sumatorio del asteroide "z" en el eje x */
				asteroides[z].sum_fX += fuerza * cos(angulo);
				/* El producto de la fuerza por el seno del angulo se agrega positvamente al sumatorio del asteroide "z" en el eje y */
				asteroides[z].sum_fY += fuerza * sin(angulo);
			}
		}
		/* Si hubo rebote contra los bordes, se posiciona el asteroide y modifica su velocidad */
		for (unsigned int j = 0; j < asteroides.size(); j++) {
			/* Calculo de las nuevas aceleraciones */
			calcularNuevaAceleracion(asteroides[j], aceleracion);

			/* Calculo de las nuevas velocidades */
			calcularNuevasVelocidades(asteroides[j], aceleracion, time_interval);

			/* Calculo de las nuevas posiciones */
			calcularNuevaPosicion(asteroides[j], time_interval);

			/* Si el asteroide rebota contra un borde, se posiciona a 5 puntos del limite del espacio
			y se modifica su velocidad */
			comprobarBordes(asteroides[j], width, height);
		}

		/* Actualizar posiciones y velocidades de cada asteroide. Los sumatorios de fuerzas se reinician a 0 */
		for (unsigned int ii = 0; ii< asteroides.size(); ii++) {
			actualizarAsteroide(asteroides[ii]);
			cout << "Velocidad del asteroide" << ii << " en iteración" << i << ":" << asteroides[ii].sig_vX << ", " << asteroides[ii].sig_vY << "\n";
		}
	}

	/* Crear el fichero de salida out.txt con un maximo de 3 decimales */
	ofstream fs2("out.txt");
	fs2 << setprecision(3) << fixed;

	/* Escribe en el fichero de salida los valores finales de cada asteroide */
	for (unsigned int i = 0; i < asteroides.size(); i++) {

		fs2 << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].vX << " " << asteroides[i].vY
				<< " " << asteroides[i].masa << endl;
	}

	/* Cierra el fichero out.txt */
	fs2.close();

	return 0;
}

///// Cuerpo de las funciones auxiliares ////
void actualizarAsteroide(asteroide ast){
	ast.pX = ast.sig_pX;
	ast.pY = ast.sig_pY;
	ast.vX = ast.sig_vX;
	ast.vY = ast.sig_vY;
	ast.sum_fX = 0;
	ast.sum_fY = 0;
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
