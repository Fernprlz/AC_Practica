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
	double pX;
	double pY;
	double masa;
	double sig_pX;
	double sig_pY;
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
void descomponerFuerzas(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo);
void calcularNuevaAceleracion(asteroide ast, double sum_fX, double sum_fY, double *aceleracion);
void calcularNuevasVelocidades(asteroide &ast, double aceleracion[2], double time_interval);
void calcularNuevaPosicion(asteroide &ast, double time_interval);
void comprobarBordes(asteroide &ast, double width, double height);
void actualizarAsteroide(asteroide *ast, double *sum_fX, double *sum_fY);
// --------------------------------------------------------------------------------------------------- //

int main(int argc, char *argv[]) {

	/* Comprueba el numero de argumentos e imprime el mensaje de error correspondiente */
	if (argc != 5 || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 || atoi(argv[3]) < 0 || atoi(argv[4]) <= 0) {

		cout << "nasteroids-par: Wrong arguments." << endl;
		cout << "Correct use:" << endl;
		cout << "nasteroids-par num_asteroides num_iteraciones num_planetas semilla" << endl;
		return -1;
	}

	int num_asteroides = atoi(argv[1]);
	int num_iteraciones = atoi(argv[2]);
	int num_planetas = atoi(argv[3]);
	int semilla = atoi(argv[4]);

	/* Definicion de constantes dadas por el enunciado del problema */
	const double gravity = 6.674e-5;
	const double time_interval = 0.1;
	const double dmin = 5.0;
	const double width = 200;
	const double height = 200;
	const double mean = 1000;
	const double sdm = 50;

	/* Random distributions */
	default_random_engine random{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width, std::numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height, std::numeric_limits<double>::max())};
	normal_distribution<double> mdist{mean, sdm};

	/* Crea un fichero de salida */
	ofstream fs("init_conf.txt");

	/* Escribe los argumentos en el fichero init_conf.txt */
	for (int ii = 0; ii < argc - 1; ii++) {

		if (ii != 3) fs << argv[ii + 1] << " ";
		else fs << argv[ii + 1] << endl;
	}

	/* Fija 3 decimales al escribir en el fichero de configuracion inicial init_conf.txt */
	fs << setprecision(3) << fixed;

	// Vector que almacena los asteroides
	vector<asteroide> asteroides(num_asteroides);
	// Vector que almacena los planetas
	vector<planeta> planetas(num_planetas);

	// Generacion de los asteroides y escritura de sus datos en el fichero init_conf.txt
	for(unsigned int ast = 0; ast < asteroides.size(); ast++) {
		asteroides[ast] = {xdist(random), ydist(random), mdist(random), 0, 0, 0, 0, 0, 0};
		fs << asteroides[ast].pX << " " << asteroides[ast].pY << " " << asteroides[ast].masa << endl;
	}


	/* Generacion de los planetas y escritura de sus datos en el fichero init_conf.txt */
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

	/* Cierra el fichero init_conf.txt */
	fs.close();

	/* Definicion e inicializacion de variables auxiliares para los calculos */
	double distancia = 0;
	double pendiente = 0;
	double angulo = 0;
	double aceleracion[2] = {0}; // 2 componentes, X e Y
	double aux_vX = 0;
	double aux_vY = 0;
	double fuerza = 0;

	/* Vectores que almacenaran en la posicion "i" el sumatorio de la fuerza del asteroide "i" en la coordenada correspondiente */
	vector<double> sum_fX(num_asteroides);
	vector<double> sum_fY(num_asteroides);

	/* Arreglos bidimensionales (modificados solo por los hilos) que almacenaran en la posicion [i][ast] la fuerza resultante
	entre el elemento "i" y "ast" en la coordenada correspondiente. Esta matriz sera compartida por todos los hilos sin
	producirse condiciones de carrera ya que ningun hilo tendra el mismo par de indices (i,ast) a la vez; a lo sumo, compartiran
	el indice "ast" pero nunca el "i" y por ende accederan a posiciones de memoria distintas */
	double **sum_total_fX = new double*[num_asteroides];
	double **sum_total_fY = new double*[num_asteroides];

	/* Matriz utilizada para marcar las colisiones entre asteroides. Un 1 en la posicion colisiones[i][ast] indica que el
	asteroide i ha colisionado con el ast. Esto para que posteriormente se recorra la matriz secuencialmente y se intercambien
	las velocidades en el orden correcto */
	double **colisiones = new double*[num_asteroides];


	/* Creacion de la segunda dimension de las matrices para cada posicion */
	for(int ast = 0; ast < num_asteroides; ast++) {

		sum_total_fX[ast] = new double[num_asteroides];
		sum_total_fY[ast] = new double[num_asteroides];
		colisiones[ast] = new double[num_asteroides];
	}

	/* Inicializacion de forma paralela de las matrices a 0 */
	#pragma omp parallel for
	for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {
		for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {
			sum_total_fX[ast_1][ast_2] = 0;
			sum_total_fY[ast_1][ast_2] = 0;
			colisiones[ast_1][ast_2] = 0;
		}
	}

	unsigned int ast_1 = 0, ast_2 = 0, pla = 0;

	for (int iteracion = 0; iteracion < num_iteraciones; iteracion++) {
		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variables definidas
		anteriormente para los calculos seran privadas para cada hilo */
		#pragma omp parallel for private(ast_2, fuerza, distancia, pendiente, angulo) /*schedule(static, num_asteroides/omp_get_num_threads())*/
		for (ast_1 = 0; ast_1 < asteroides.size(); ast_1++) {

			/* Se evalua la interaccion entre el asteroide "ast" y los demas que tengan índice mayor (a partir de ast+1) */
			for (ast_2 = ast_1 + 1; ast_2 < asteroides.size(); ast_2++) {

				/* Calcula la distancia entre el asteroide "ast" y el asteroide "k" */
				distancia = calcularDistanciaAsteroide(asteroides[ast_1], asteroides[ast_2]);

				//TODO: ESTO ES LO QUE DECIA EL PIBE QUE EL NASTEROIDS NO CUMPLE, QUE EL ASTEROIDE 0 NO LO USA PANADA
				/* Si la distancia es menor o igual a 2 se marca con un 1 que el asteroide ast colisiono con el asteroide k */
				if (distancia <= dmin) {
					if (ast_1 > 0) {
						colisiones[ast_1][ast_2] = 1;
					}
				} else {
					/* Si las distancia es mayor a la distancia minima (5) se toma encuenta la fuerza de atraccion entre ellos */

					/* Calcula la pendiente entre el asteroide "ast" y el asteroide "k" */
					pendiente =  calcularPendienteAsteroide(asteroides[ast_1], asteroides[ast_2]);

					/* El angulo entre asteroides sera la ArcTan(pendiente) */
					angulo = atan(pendiente);

					/* Calculo de la fuerza resultante entre ambos asteroides */
					fuerza = calcularFuerzaAsteroide(asteroides[ast_1], asteroides[ast_2], gravity, distancia);

					/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje x */				
					descomponerFuerzas(&sum_total_fX[ast_1][ast_2], &sum_total_fX[ast_2][ast_1], fuerza, angulo);
					/* La fuerza resultante por el seno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje y */
					descomponerFuerzas(&sum_total_fY[ast_1][ast_2], &sum_total_fY[ast_2][ast_1], fuerza, angulo);
				}
			}
		}

		/* Se copia paralelamente a la posicion "i" de los vectores sum_fX y sum_fY el sumatorio de toda la fila
		"i" de la matriz de la coordenada correspondiente. Dicho sumatorio se corresponde con las fuerza totales del asteroide
		"i" para con los demas elementos */
		#pragma omp parallel for
		for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {

			for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {

				sum_fX[ast_1] += sum_total_fX[ast_1][ast_2];
				sum_fY[ast_1] += sum_total_fY[ast_1][ast_2];
			}
		}

		/* Bucle anidado que se ejecuta de forma secuencial ya que por estar accediendo a las velocidades de los asteroides
		no es posible paralelizarlo. Se recorre en orden (fila a fila) la matriz "colisiones" para gestionar los choques en
		orden correcto. Si colisiones[i][ast] == 1 se intercambia la velocidad del asteroide i con el del ast, y asi sucesivamente */
		for (int ast_1 = 0; ast_1 < num_asteroides; ast_1++) {

			for (int ast_2 = ast_1 + 1; ast_2 < num_asteroides; ast_2++) {

				if (colisiones[ast_1][ast_2] == 1) {
					aux_vX = asteroides[ast_1].vX;
					aux_vY = asteroides[ast_1].vY;
					asteroides[ast_1].vX = asteroides[ast_2].vX;
					asteroides[ast_1].vY = asteroides[ast_2].vY;
					asteroides[ast_2].vX = aux_vX;
					asteroides[ast_2].vY = aux_vY;
				}
			}
		}

		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variblaes definidas
		anteriormente para los calculos seran privadas para cada hilo. A diferencia del doble bucle anterior, en este se
		evalua la fuerza entre los asteroides y planetas, almacenando la fuerza resultante de cada interaccion en los
		vectores unidimensionales rellenados recientemente. Ningun hilo accedera a la misma posicion de un vector ya que a lo
		sumo podran coincidir en estar evaluando el mismo planeta (indice "pla") pero nunca el mismo asteroide (con lo cual
		el indice "ast" sera distinto) */
		#pragma omp parallel for private(pla, fuerza, distancia, pendiente, angulo)
		for (unsigned int ast = 0; ast < asteroides.size(); ast++) {

			/* Para un mismo asteroide "ast" se itera sobre todos los planetas "y" evaluando la interaccion para cada par de
			elementos */
			for (pla = 0; pla < planetas.size(); pla++) {

				/* Calcula la distancia entre el asteroide "ast" y el planeta "y" */
				distancia = calcularDistanciaPlaneta(planetas[pla], asteroides[ast]);
				/* Calculo de la pendiente entre elementos */
			  pendiente = calcularPendientePlaneta(planetas[pla], asteroides[ast]);

				/* EL angulo entre asteroide y planeta sera el ArcTan(pendiente) */
				angulo = atan(pendiente);

				// TODO: FUERZA UMBRAL EN LA QUE TRUNCAR OK
				/* Calculo de la fuerza en el resultante entre ambos elementos */
				fuerza = calcularFuerzaPlaneta(planetas[pla], asteroides[ast], gravity, distancia);

				/* El producto de la fuerza por el coseno del angulo se agrega positvamente al sumatorio del asteroide "ast" en el eje x */
				sum_fX[ast] += fuerza * cos(angulo);

				/* El producto de la fuerza por el seno del angulo se agrega positvamente al sumatorio del asteroide "ast" en el eje y */
				sum_fY[ast] += fuerza * sin(angulo);
			}
		}

		/* Bucle en el que se realizan todos los calculos para obtener las nuevas velocidades y posiciones de cada asteroide */
		/* Cada hilo tendra su copia de la variable aceleracionx y aceleraciony. */
		#pragma omp parallel for private(aceleracion)
		for (unsigned int ast = 0; ast < asteroides.size(); ast++) {

			/* Calculo de las nuevas aceleraciones ->>>>> HAY QUE METER SUMFX TB*/
			calcularNuevaAceleracion(asteroides[ast], sum_fX[ast], sum_fY[ast], aceleracion);

			/* Calculo de las nuevas velocidades*/
			calcularNuevasVelocidades(asteroides[ast], aceleracion, time_interval);

			/* Calculo de las nuevas posiciones*/
			calcularNuevaPosicion(asteroides[ast], time_interval);

			/* Se hubo rebote contra los bordes, se reposiciona el asteroide y modifica su velocidad */
			comprobarBordes(asteroides[ast], width, height);
		}

		/* Para cada asteroide se actualizan paralelamente las nuevas posiciones y velocidades. Tambien, todos los sumatorios
		se actualizan a 0 */
		#pragma omp parallel for
		for (unsigned int ast_1 = 0; ast_1 < asteroides.size(); ast_1++) {

			actualizarAsteroide(&asteroides[ast_1], &sum_fX[ast_1], &sum_fY[ast_1]);

			for (int ast_2 = 0; ast_2 < num_asteroides; ast_2++) {

				sum_total_fX[ast_1][ast_2] = 0;
				sum_total_fY[ast_1][ast_2] = 0;
				colisiones[ast_1][ast_2] = 0;
			}

		}
		for(unsigned int ii = 0; ii < asteroides.size(); ii++){
			cout << "Velocidad del asteroide [" << ii << "]: " << asteroides[ii].sig_vX << ", " << asteroides[ii].sig_vY << "\n";
		}
	}

	/* Crea el fichero de salida out.txt y fija 3 decimales al escribir en el */
	ofstream fs2("out.txt");
	fs2 << setprecision(3) << fixed;

	/* Escribe en el fichero de salida los valores finales de cada asteroide */
	for (unsigned int i = 0; i < asteroides.size(); i++) {

		fs2 << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].vX << " " << asteroides[i].vY
		<< " " << asteroides[i].masa << endl;
	}

	/* Cierra el fichero out.txt */
	fs2.close();

	/* Libera todas la memoria reservada para las matrices que almacenan las fuerzas */
	#pragma omp parallel for
	for(int i = 0 ; i < num_asteroides; i++) {

		delete[] sum_total_fX[i];
		delete[] sum_total_fY[i];
		delete[] colisiones[i];
	}

	delete[] sum_total_fX;
	delete[] sum_total_fY;
	delete[] colisiones;

	return 0;
}




/////////////////////////////// Cuerpo de las funciones auxiliares ///////////////////////////////////
double calcularDistanciaAsteroide(asteroide cuerpo1, asteroide cuerpo2){
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
double calcularFuerzaAsteroide(asteroide cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
	/* Si la fuerza resultante es mayor a 100, se trunca a este valor */
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
double calcularFuerzaPlaneta(planeta cuerpo1, asteroide cuerpo2, double gravity, double distancia){
	double fuerza = (gravity * cuerpo1.masa * cuerpo2.masa) / (pow(distancia, 2.0));
		/* Si la fuerza resultante es mayor a 100, se trunca a este valor */
		if (fuerza > 100){
			fuerza = 100;
		}
		return fuerza;
}

void descomponerFuerzas(double *sum_total_f1, double *sum_total_f2, double fuerza, double angulo){
	/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide ast y se resta al asteroide k en el eje x */
	*sum_total_f1 += fuerza * cos(angulo);
	*sum_total_f2 -= fuerza * cos(angulo);
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
