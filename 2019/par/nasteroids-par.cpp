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

struct asteroide {
	double posx;
	double posy;
	double masa;
	double nueva_posx;
	double nueva_posy;
	double vx;
	double vy;
	double nueva_vx;
	double nueva_vy;
	int colision;
	int id_colision_desde_inicio;
};

struct planeta {
	double posx;
	double posy;
	double masa;
};

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
	const double gravedad = 6.674e-5;
	const double intervalo_tiempo = 0.1;
	const double distancia_min = 5.0;
	const double ancho = 200;
	const double alto = 200;
	const double media = 1000;
	const double desviacion = 50;

	/* Random distributions */
	default_random_engine re{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(ancho, std::numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(alto, std::numeric_limits<double>::max())};
	normal_distribution<double> mdist{media, desviacion};

	/* Crea un fichero de salida */
   	ofstream fs("init_conf.txt");

   	/* Escribe los argumentos en el fichero init_conf.txt */
   	for (int i = 0; i < argc - 1; i++) {

   		if (i != 3) fs << argv[i+1] << " ";
   		else fs << argv[i+1] << endl;
   	}

   	/* Fija 3 decimales al escribir en el fichero de configuracion inicial init_conf.txt */
   	fs << setprecision(3) << fixed;

	vector<asteroide> asteroides(num_asteroides); /* Vector que almacena los asteroides */

	/* Generacion de los asteroides y escritura de sus datos en el fichero init_conf.txt */
	for(unsigned int i = 0; i < asteroides.size(); i++) {

		asteroides[i] = {xdist(re), ydist(re), mdist(re), 0, 0, 0, 0, 0, 0, 0, 0};
		fs << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].masa << endl;
	}

	vector<planeta> planetas(num_planetas); /* Vector que almacena los planetas */

	/* Generacion de los planetas y escritura de sus datos en el fichero init_conf.txt */
	for(unsigned int i = 0; i < planetas.size(); i++) {

		if (i % 4 == 0) {
			planetas[i] = {0.0, ydist(re), mdist(re) * 10};
		}
		else if (i % 4 == 1) {
			planetas[i] = {xdist(re), 0.0, mdist(re) * 10};
		}
		else if (i % 4 == 2) {
			planetas[i] = {ancho, ydist(re), mdist(re) * 10};
		}
		else if (i % 4 == 3) {
			planetas[i] = {xdist(re), alto, mdist(re) * 10};
		}
		fs << planetas[i].posx << " " << planetas[i].posy << " " << planetas[i].masa << endl;
	}

	/* Cierra el fichero init_conf.txt */
	fs.close();

	/* Definicion e inicializacion de variables auxiliares para los calculos */
	double distancia_ab = 0;
	double pendiente_ab = 0;
	double angulo_ab = 0;
	double aceleracionx = 0;
	double aceleraciony = 0;
	double vxaux = 0;
	double vyaux = 0;
	double fuerza = 0;

	/* Vectores que almacenaran en la posicion "i" el sumatorio de la fuerza del asteroide "i" en la coordenada correspondiente */
	vector<double> sumatorios_fx(num_asteroides);
	vector<double> sumatorios_fy(num_asteroides);

	/* Arreglos bidimensionales (modificados solo por los hilos) que almacenaran en la posicion [i][j] la fuerza resultante
	entre el elemento "i" y "j" en la coordenada correspondiente. Esta matriz sera compartida por todos los hilos sin
	producirse condiciones de carrera ya que ningun hilo tendra el mismo par de indices (i,j) a la vez; a lo sumo, compartiran
	el indice "j" pero nunca el "i" y por ende accederan a posiciones de memoria distintas */
	double **sumatorioX = new double*[num_asteroides];
	double **sumatorioY = new double*[num_asteroides];

	/* Matriz utilizada para marcar las colisiones entre asteroides. Un 1 en la posicion colisiones[i][j] indica que el
	asteroide i ha colisionado con el j. Esto para que posteriormente se recorra la matriz secuencialmente y se intercambien
	las velocidades en el orden correcto */
	double **colisiones = new double*[num_asteroides];


	/* Creacion de la segunda dimension de las matrices para cada posicion */
	for(int i = 0; i < num_asteroides; i++) {

		 sumatorioX[i] = new double[num_asteroides];
		 sumatorioY[i] = new double[num_asteroides];
		 colisiones[i] = new double[num_asteroides];
	}

	/* Inicializacion de forma paralela de las matrices a 0 */
	#pragma omp parallel for
	for (int i = 0; i < num_asteroides; i++) {

		for (int j = 0; j < num_asteroides; j++) {

			sumatorioX[i][j] = 0;
			sumatorioY[i][j] = 0;
			colisiones[i][j] = 0;
		}
	}

	unsigned int k = 0, j = 0, y = 0;

	for (int i = 0; i < num_iteraciones; i++) {

		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variblaes definidas
		anteriormente para los calculos seran privadas para cada hilo */
		#pragma omp parallel for private(k, fuerza, distancia_ab, pendiente_ab, angulo_ab) /*schedule(static, num_asteroides/omp_get_num_threads())*/
		for (j = 0; j < asteroides.size(); j++) {

			/* Se evalua la interaccion entre el asteroide "j" y los demas que tengan índice mayor (a partir de j+1) */
			for (k = j+1; k < asteroides.size(); k++) {

				// TODO: DISTANCIA OK
				/* Calcula la distancia entre el asteroide "j" y el asteroide "k" */
				distancia_ab = sqrt(pow(asteroides[j].posx - asteroides[k].posx, 2.0) +
									pow(asteroides[j].posy - asteroides[k].posy, 2.0));

				/* Si las distancia es mayor a la distancia minima (5) se toma encuenta la fuerza de atraccion entre ellos */
				if (distancia_ab > distancia_min) {

					// TODO: PENDIENTE OK
					/* Calcula la pendiente entre el asteroide "j" y el asteroide "k" */
					pendiente_ab = (asteroides[j].posy - asteroides[k].posy) / (asteroides[j].posx - asteroides[k].posx);

					// TODO: CORRECCION DE LA PENDIENTE
					/* Si la pendiente es mayor a 1 o menor a -1 se modifica el valor de la misma */
					if (pendiente_ab > 1){
							pendiente_ab = 1;
					} else if (pendiente_ab < -1) {
							pendiente_ab = -1;
					}

					/* El angulo entre asteroides sera la ArcTan(pendiente) */
					angulo_ab = atan(pendiente_ab);

					/* Calculo de la fuerza resultante entre ambos asteroides */
					fuerza = ((gravedad * asteroides[j].masa * asteroides[k].masa) / (pow(distancia_ab, 2.0)));

					/* Si la fuerza resultante es mayor a 200, se trunca a este valor */
					if (fuerza > 200) fuerza = 200;

					/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje x */
					sumatorioX[j][k] += fuerza * cos(angulo_ab); /* Al asteroide "j" se le suma la fuerza calculada en el eje x */
					sumatorioX[k][j] -= fuerza * cos(angulo_ab); /* Al asteroide "k" se le resta la fuerza calculada en el eje x */

					/* La fuerza resultante por el seno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje y */
					sumatorioY[j][k] += fuerza * sin(angulo_ab); /* Al asteroide "j" se le suma la fuerza calculada en el eje y */
					sumatorioY[k][j] -= fuerza * sin(angulo_ab); /* Al asteroide "k" se le resta la fuerza calculada en el eje y */
				}
				/* Si la distancia es menor o igual a 2 se marca con un 1 que el asteroide j colisiono con el asteroide k */
				else {

					if (j > 0) {

						colisiones[j][k] = 1;
					}
				}
			}
		}

		/* Se copia paralelamente a la posicion "i" de los vectores sumatorios_fx y sumatorios_fy el sumatorio de toda la fila
		"i" de la matriz de la coordenada correspondiente. Dicho sumatorio se corresponde con las fuerza totales del asteroide
		"i" para con los demas elementos */
		#pragma omp parallel for
		for (int i = 0; i < num_asteroides; i++) {

			for (int j = 0; j < num_asteroides; j++) {

				sumatorios_fx[i] += sumatorioX[i][j];
				sumatorios_fy[i] += sumatorioY[i][j];
			}
		}

		/* Bucle anidado que se ejecuta de forma secuencial ya que por estar accediendo a las velocidades de los asteroides
		no es posible paralelizarlo. Se recorre en orden (fila a fila) la matri "colisiones" para gestionar los choques en
		orden correcto. Si colisiones[i][j] == 1 se intercambia la velocidad del asteroide i con el del j, y asi sucesivamente */
		for (int i = 0; i < num_asteroides; i++) {

			for (int j = i+1; j < num_asteroides; j++) {

				if (colisiones[i][j] == 1) {
					vxaux = asteroides[i].vx;
					vyaux = asteroides[i].vy;
					asteroides[i].vx = asteroides[j].vx;
					asteroides[i].vy = asteroides[j].vy;
					asteroides[j].vx = vxaux;
					asteroides[j].vy = vyaux;
				}
			}
		}

		/* Cada hilo ejecutara cierto numero de iteraciones del bucle anidado a continuacion. Las variblaes definidas
		anteriormente para los calculos seran privadas para cada hilo. A diferencia del doble bucle anterior, en este se
		evalua la fuerza entre los asteroides y planetas, almacenando la fuerza resultante de cada interaccion en los
		vectores unidimensionales rellenados recientemente. Ningun hilo accedera a la misma posicion de un vector ya que a lo
		sumo podran coincidir en estar evaluando el mismo planeta (indice "y") pero nunca el mismo asteroide (con lo cual
		el indice "j" sera distinto) */
		#pragma omp parallel for private(y, fuerza, distancia_ab, pendiente_ab, angulo_ab)
		for (unsigned int j = 0; j < asteroides.size(); j++) {

			/* Para un mismo asteroide "j" se itera sobre todos los planetas "y" evaluando la interaccion para cada par de
			elementos */
			for (y = 0; y < planetas.size(); y++) {

				/* Calcula la distancia entre el asteroide "j" y el planeta "y" */
				distancia_ab = sqrt(pow(planetas[y].posx - asteroides[j].posx, 2.0) +
									pow(planetas[y].posy - asteroides[j].posy, 2.0));

				/* Calculo de la pendiente entre elementos */
				pendiente_ab = (planetas[y].posy - asteroides[j].posy) / (planetas[y].posx - asteroides[j].posx);

				/* Si la pendiente es mayor a 1 o menor a -1 se modifica el valor de la misma */
				if (pendiente_ab > 1 || pendiente_ab < -1) pendiente_ab -= int(pendiente_ab / 1);

				/* EL angulo entre asteroide y planeta sera el ArcTan(pendiente) */
				angulo_ab = atan(pendiente_ab);

				// TODO: FUERZA UMBRAL EN LA QUE TRUNCAR OK
				/* Calculo de la fuerza en el resultante entre ambos elementos */
				fuerza = ((gravedad * planetas[y].masa * asteroides[j].masa) / (pow(distancia_ab, 2.0)));
				/* Si la fuerza resultante es mayor a 100 se trunca a este valor */
				if (fuerza > 100) fuerza = 100;

				/* El producto de la fuerza por el coseno del angulo se agrega positvamente al sumatorio del asteroide "j" en el eje x */
				sumatorios_fx[j] += fuerza * cos(angulo_ab);

				/* El producto de la fuerza por el seno del angulo se agrega positvamente al sumatorio del asteroide "j" en el eje y */
				sumatorios_fy[j] += fuerza * sin(angulo_ab);
			}
		}

		/* Bucle en el que se realizan todos los calculos para obtener las nuevas velocidades y posiciones de cada asteroide */
		/* Cada hilo tendra su copia de la variable aceleracionx y aceleraciony. */
		#pragma omp parallel for private(aceleracionx, aceleraciony)
		for (unsigned int j = 0; j < asteroides.size(); j++) {

			aceleracionx = 0;
			aceleraciony = 0;

			/* Calculo de las nuevas aceleraciones */
			aceleracionx = sumatorios_fx[j] / asteroides[j].masa;
			aceleraciony = sumatorios_fy[j] / asteroides[j].masa;

			/* Calculo de las nuevas velocidades */
			asteroides[j].nueva_vx = asteroides[j].vx + (aceleracionx * intervalo_tiempo);
			asteroides[j].nueva_vy = asteroides[j].vy + (aceleraciony * intervalo_tiempo);

			/* Calculo de las nuevas posiciones */
			asteroides[j].nueva_posx = asteroides[j].posx + (asteroides[j].nueva_vx * intervalo_tiempo);
			asteroides[j].nueva_posy = asteroides[j].posy + (asteroides[j].nueva_vy * intervalo_tiempo);

			// TODO: REBOTE CON EL BORDE OK
			/* Se hubo rebote contra los bordes, se reposiciona el asteroide y modifica su velocidad */
			if (asteroides[j].nueva_posx <= 0) {
				asteroides[j].nueva_posx = 5;
				asteroides[j].nueva_vx *= (-1);
			}
			if (asteroides[j].nueva_posx >= ancho) {
				asteroides[j].nueva_posx = ancho - 5;
				asteroides[j].nueva_vx *= (-1);
			}
			if (asteroides[j].nueva_posy <= 0) {
				asteroides[j].nueva_posy = 5;
				asteroides[j].nueva_vy *= (-1);
			}
			if (asteroides[j].nueva_posy >= alto) {
				asteroides[j].nueva_posy = alto - 5;
				asteroides[j].nueva_vy *= (-1);
			}
		}

		/* Para cada asteroide se actualizan paralelamente las nuevas posiciones y velocidades. Tambien, todos los sumatorios
		se actualizan a 0 */
		#pragma omp parallel for
		for (unsigned int i = 0; i< asteroides.size(); i++) {

			asteroides[i].posx = asteroides[i].nueva_posx;
			asteroides[i].posy = asteroides[i].nueva_posy;
			asteroides[i].vx = asteroides[i].nueva_vx;
			asteroides[i].vy = asteroides[i].nueva_vy;
			sumatorios_fx[i] = 0;
			sumatorios_fy[i] = 0;

			for (int j = 0; j < num_asteroides; j++) {

				sumatorioX[i][j] = 0;
				sumatorioY[i][j] = 0;
				colisiones[i][j] = 0;
			}
		}
	}

	/* Crea el fichero de salida out.txt y fija 3 decimales al escribir en el */
	ofstream fs2("out.txt");
	fs2 << setprecision(3) << fixed;

	/* Escribe en el fichero de salida los valores finales de cada asteroide */
	for (unsigned int i = 0; i < asteroides.size(); i++) {

		fs2 << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].vx << " " << asteroides[i].vy
				<< " " << asteroides[i].masa << endl;
	}

	/* Cierra el fichero out.txt */
	fs2.close();

	/* Libera todas la memoria reservada para las matrices que almacenan las fuerzas */
	#pragma omp parallel for
	for(int i = 0 ; i < num_asteroides; i++) {

		delete[] sumatorioX[i];
		delete[] sumatorioY[i];
		delete[] colisiones[i];
	}

	delete[] sumatorioX;
	delete[] sumatorioY;
	delete[] colisiones;

	return 0;
}
