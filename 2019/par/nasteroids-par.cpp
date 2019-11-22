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
	double sig_vy;
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
	for (int i = 0; i < argc - 1; i++) {

		if (i != 3) fs << argv[i+1] << " ";
		else fs << argv[i+1] << endl;
	}

	/* Fija 3 decimales al escribir en el fichero de configuracion inicial init_conf.txt */
	fs << setprecision(3) << fixed;

	// Vector que almacena los asteroides
	vector<asteroide> asteroides(num_asteroides);
	// Vector que almacena los planetas
	vector<planeta> planetas(num_planetas);

	// Generacion de los asteroides y escritura de sus datos en el fichero init_conf.txt
	for(unsigned int i = 0; i < asteroides.size(); i++) {
		asteroides[i] = {xdist(random), ydist(random), mdist(random), 0, 0, 0, 0, 0, 0, 0, 0};
		fs << asteroides[i].pX << " " << asteroides[i].pY << " " << asteroides[i].masa << endl;
	}



	/* Generacion de los planetas y escritura de sus datos en el fichero init_conf.txt */
	for(unsigned int i = 0; i < planetas.size(); i++) {

		if (i % 4 == 0) {
			planetas[i] = {0.0, ydist(random), mdist(random) * 10};
		}
		else if (i % 4 == 1) {
			planetas[i] = {xdist(random), 0.0, mdist(random) * 10};
		}
		else if (i % 4 == 2) {
			planetas[i] = {width, ydist(random), mdist(random) * 10};
		}
		else if (i % 4 == 3) {
			planetas[i] = {xdist(random), height, mdist(random) * 10};
		}
		fs << planetas[i].pX << " " << planetas[i].pY << " " << planetas[i].masa << endl;
	}

	/* Cierra el fichero init_conf.txt */
	fs.close();

	/* Definicion e inicializacion de variables auxiliares para los calculos */
	double distancia = 0;
	double pendiente = 0;
	double angulo = 0;
	double aceleracion[2] = {0}; // 2 componentes, X e Y
	double v_Xaux = 0;
	double vY_aux = 0;
	double fuerza = 0;

	/* Vectores que almacenaran en la posicion "i" el sumatorio de la fuerza del asteroide "i" en la coordenada correspondiente */
	vector<double> sum_fX(num_asteroides);
	vector<double> sum_fY(num_asteroides);

	/* Arreglos bidimensionales (modificados solo por los hilos) que almacenaran en la posicion [i][j] la fuerza resultante
	entre el elemento "i" y "j" en la coordenada correspondiente. Esta matriz sera compartida por todos los hilos sin
	producirse condiciones de carrera ya que ningun hilo tendra el mismo par de indices (i,j) a la vez; a lo sumo, compartiran
	el indice "j" pero nunca el "i" y por ende accederan a posiciones de memoria distintas */
	double **sum_total_fX = new double*[num_asteroides];
	double **sum_total_fY = new double*[num_asteroides];

	/* Matriz utilizada para marcar las colisiones entre asteroides. Un 1 en la posicion colisiones[i][j] indica que el
	asteroide i ha colisionado con el j. Esto para que posteriormente se recorra la matriz secuencialmente y se intercambien
	las velocidades en el orden correcto */
	double **colisiones = new double*[num_asteroides];


	/* Creacion de la segunda dimension de las matrices para cada posicion */
	for(int i = 0; i < num_asteroides; i++) {

		sum_total_fX[i] = new double[num_asteroides];
		sum_total_fY[i] = new double[num_asteroides];
		colisiones[i] = new double[num_asteroides];
	}

	/* Inicializacion de forma paralela de las matrices a 0 */
	#pragma omp parallel for
	for (int i = 0; i < num_asteroides; i++) {
		for (int j = 0; j < num_asteroides; j++) {
			sum_total_fX[i][j] = 0;
			sum_total_fY[i][j] = 0;
			colisiones[i][j] = 0;
		}
	}


	unsigned int k = 0, j = 0, y = 0;

	for (int i = 0; i < num_iteraciones; i++) {

		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variblaes definidas
		anteriormente para los calculos seran privadas para cada hilo */
		#pragma omp parallel for private(k, fuerza, distancia, pendiente, angulo) /*schedule(static, num_asteroides/omp_get_num_threads())*/
		for (j = 0; j < asteroides.size(); j++) {

			/* Se evalua la interaccion entre el asteroide "j" y los demas que tengan índice mayor (a partir de j+1) */
			for (k = j+1; k < asteroides.size(); k++) {

				/* Calcula la distancia entre el asteroide "j" y el asteroide "k" */
				distancia = calcularDistanciaAsteroide(asteroides[j], asteroides[k]);

				/* Si la distancia es menor o igual a 2 se marca con un 1 que el asteroide j colisiono con el asteroide k */
				if (distancia <= dmin) {
					if (j > 0) {
						colisiones[j][k] = 1;
					}
				} else {
					/* Si las distancia es mayor a la distancia minima (5) se toma encuenta la fuerza de atraccion entre ellos */

					/* Calcula la pendiente entre el asteroide "j" y el asteroide "k" */
					pendiente =  calcularPendienteAsteroide(asteroides[j], asteroides[k]);

					/* El angulo entre asteroides sera la ArcTan(pendiente) */
					angulo = atan(pendiente);

					/* Calculo de la fuerza resultante entre ambos asteroides */
					fuerza = calcularFuerzaAsteroide(asteroides[j], asteroides[k], gravity, distancia);

					// TODO: Cambiar funcion de descomponer fuerzas.
					/////////////////////
					/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje x */
					sum_total_fX[j][k] += fuerza * cos(angulo); /* Al asteroide "j" se le suma la fuerza calculada en el eje x */
					sum_total_fX[k][j] -= fuerza * cos(angulo); /* Al asteroide "k" se le resta la fuerza calculada en el eje x */

					/* La fuerza resultante por el seno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje y */
					sum_total_fY[j][k] += fuerza * sin(angulo); /* Al asteroide "j" se le suma la fuerza calculada en el eje y */
					sum_total_fY[k][j] -= fuerza * sin(angulo); /* Al asteroide "k" se le resta la fuerza calculada en el eje y */
					//////////////////////
				}
			}
		}

		/* Se copia paralelamente a la posicion "i" de los vectores sum_fX y sum_fY el sumatorio de toda la fila
		"i" de la matriz de la coordenada correspondiente. Dicho sumatorio se corresponde con las fuerza totales del asteroide
		"i" para con los demas elementos */
		#pragma omp parallel for
		for (int i = 0; i < num_asteroides; i++) {

			for (int j = 0; j < num_asteroides; j++) {

				sum_fX[i] += sum_total_fX[i][j];
				sum_fY[i] += sum_total_fY[i][j];
			}
		}

		/* Bucle anidado que se ejecuta de forma secuencial ya que por estar accediendo a las velocidades de los asteroides
		no es posible paralelizarlo. Se recorre en orden (fila a fila) la matriz "colisiones" para gestionar los choques en
		orden correcto. Si colisiones[i][j] == 1 se intercambia la velocidad del asteroide i con el del j, y asi sucesivamente */
		for (int i = 0; i < num_asteroides; i++) {

			for (int j = i+1; j < num_asteroides; j++) {

				if (colisiones[i][j] == 1) {
					vXaux = asteroides[i].vX;
					vYaux = asteroides[i].vY;
					asteroides[i].vX = asteroides[j].vX;
					asteroides[i].vY = asteroides[j].vY;
					asteroides[j].vX = vXaux;
					asteroides[j].vY = vYaux;
				}
			}
		}

		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variblaes definidas
		anteriormente para los calculos seran privadas para cada hilo. A diferencia del doble bucle anterior, en este se
		evalua la fuerza entre los asteroides y planetas, almacenando la fuerza resultante de cada interaccion en los
		vectores unidimensionales rellenados recientemente. Ningun hilo accedera a la misma posicion de un vector ya que a lo
		sumo podran coincidir en estar evaluando el mismo planeta (indice "y") pero nunca el mismo asteroide (con lo cual
		el indice "j" sera distinto) */
		#pragma omp parallel for private(y, fuerza, distancia, pendiente, angulo)
		for (unsigned int j = 0; j < asteroides.size(); j++) {

			/* Para un mismo asteroide "j" se itera sobre todos los planetas "y" evaluando la interaccion para cada par de
			elementos */
			for (y = 0; y < planetas.size(); y++) {

				/* Calcula la distancia entre el asteroide "j" y el planeta "y" */
				distancia = calcularDistanciaPlaneta(planetas[q], asteroides[z]);

				/* Calculo de la pendiente entre elementos */
			  pendiente = calcularPendientePlaneta(planetas[q], asteroides[z]);

				/* EL angulo entre asteroide y planeta sera el ArcTan(pendiente) */
				angulo = atan(pendiente);

				// TODO: FUERZA UMBRAL EN LA QUE TRUNCAR OK
				/* Calculo de la fuerza en el resultante entre ambos elementos */
				fuerza = calcularFuerzaPlaneta(planetas[q], asteroides[z], gravity, distancia);

				/* El producto de la fuerza por el coseno del angulo se agrega positvamente al sumatorio del asteroide "j" en el eje x */
				sum_fX[j] += fuerza * cos(angulo);

				/* El producto de la fuerza por el seno del angulo se agrega positvamente al sumatorio del asteroide "j" en el eje y */
				sum_fY[j] += fuerza * sin(angulo);
			}
		}

		/* Bucle en el que se realizan todos los calculos para obtener las nuevas velocidades y posiciones de cada asteroide */
		/* Cada hilo tendra su copia de la variable aceleracionx y aceleraciony. */
		#pragma omp parallel for private(aceleracion)
		for (unsigned int j = 0; j < asteroides.size(); j++) {
			/* TODO: En el secuencial uso un array de 2  posiciones, a ver si no jode el pragma
			aceleracionx = 0;
			aceleraciony = 0;
			*/
			/* Calculo de las nuevas aceleraciones ->>>>> HAY QUE METER SUMFX TB
			//aceleracionx = sum_fX[j] / asteroides[j].masa;
			//aceleraciony = sum_fY[j] / asteroides[j].masa;*/
			calcularNuevaAceleracion(asteroides[j], sum_fX[j], sum_fY[j], aceleracion);

			/* Calculo de las nuevas velocidades
			asteroides[j].sig_vX = asteroides[j].vX + (aceleracionx * time_interval);
			asteroides[j].sig_vy = asteroides[j].vY + (aceleraciony * time_interval);*/
			calcularNuevasVelocidades(asteroides[j], aceleracion, time_interval);

			/* Calculo de las nuevas posiciones
			asteroides[j].sig_pX = asteroides[j].pX + (asteroides[j].sig_vX * time_interval);
			asteroides[j].sig_pY = asteroides[j].pY + (asteroides[j].sig_vy * time_interval);*/
			calcularNuevaPosicion(asteroides[j], time_interval);

			// TODO: REBOTE CON EL BORDE OK
			/* Se hubo rebote contra los bordes, se reposiciona el asteroide y modifica su velocidad */
			comprobarBordes(asteroides[j], width, height);
		}

		/* Para cada asteroide se actualizan paralelamente las nuevas posiciones y velocidades. Tambien, todos los sumatorios
		se actualizan a 0 */
		#pragma omp parallel for
		for (unsigned int i = 0; i< asteroides.size(); i++) {

		/*	asteroides[i].pX = asteroides[i].sig_pX;
			asteroides[i].pY = asteroides[i].sig_pY;
			asteroides[i].vX = asteroides[i].sig_vX;
			asteroides[i].vY = asteroides[i].sig_vy;
			sum_fX[i] = 0; 				TODO: OJO CON ESTOS DOS CAMPOS, EN EL SECUENCIAL SON PARTE DEL ASTEROIDE, AQUI ES UNA VARIABLE APARTE
			sum_fY[i] = 0;*/
			actualizarAsteroide(&asteroides[ii]);

			for (int j = 0; j < num_asteroides; j++) {

				sum_total_fX[i][j] = 0;
				sum_total_fY[i][j] = 0;
				colisiones[i][j] = 0;
			}
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
