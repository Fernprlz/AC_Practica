#include <iostream>
#include <random>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

struct asteroide {
	double posx;
	double posy;
	double masa;
	double nueva_posx;
	double nueva_posy;
	double sum_fx;
	double sum_fy;
	double vx;
	double vy;
	double nueva_vx;
	double nueva_vy;
};

struct planeta {
	double posx;
	double posy;
	double masa;
};

int main(int argc, char *argv[]) {

	/* Comprueba el numero de argumentos e imprime el mensaje de error correspondiente */
	if (argc != 5 || atoi(argv[1]) < 0 || atoi(argv[2]) < 0 || atoi(argv[3]) < 0 || atoi(argv[4]) <= 0) {

		cout << "nasteroids-seq: Wrong arguments." << endl;
		cout << "Correct use:" << endl;
		cout << "nasteroids-seq num_asteroides num_iteraciones num_planetas semilla" << endl;
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

	// Random distributions
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

   	/* Fija 3 decimales al escribir en los ficheros */
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
			planetas[i] = {200.0, ydist(re), mdist(re) * 10};
		}
		else if (i % 4 == 3) {
			planetas[i] = {xdist(re), 200.0, mdist(re) * 10};
		}
		fs << planetas[i].posx << " " << planetas[i].posy << " " << planetas[i].masa << endl;
	}

	/* Cierra el fichero init_cong.txt */
	fs.close();

	/* Definicial e inicializacion de variables auxiliares para los calculos */
	double distancia_ab = 0;
	double pendiente_ab = 0;
	double angulo_ab = 0;
	double aceleracionx = 0;
	double aceleraciony = 0;
	double vxaux = 0;
	double vyaux = 0;
	double fuerza = 0;

	for (int i = 0; i < num_iteraciones; i++) {

		for (unsigned int j = 0; j < asteroides.size(); j++) {

			/* Se evalua la interaccion entre el asteroide "j" y los demas que tengan índice mayor (a partir de j+1) */
			for (unsigned int k = j+1; k < asteroides.size(); k++) {

				/* Calcula la distancia entre el asteroide "j" y el asteroide "k" */
				distancia_ab = sqrt(pow(asteroides[j].posx - asteroides[k].posx, 2.0) +
									pow(asteroides[j].posy - asteroides[k].posy, 2.0));

				/* Si la distancia entre ambos asteroides es menor a 2 se intercambian las velocidades de la iteracion anterior */
				if (distancia_ab <= distancia_min) {

					if (j > 0) {

						vxaux = asteroides[j].vx;
						vyaux = asteroides[j].vy;
						asteroides[j].vx = asteroides[k].vx;
						asteroides[j].vy = asteroides[k].vy;
						asteroides[k].vx = vxaux;
						asteroides[k].vy = vyaux;
					}
				}

				/* Si las distancia es mayor a la distancia minima (2) se toma en cuenta la fuerza de atraccion entre ellos */
				else {

					/* Calcula la pendiente entre el asteroide "j" y el asteroide "k" */
					pendiente_ab = (asteroides[j].posy - asteroides[k].posy) / (asteroides[j].posx - asteroides[k].posx);

					/* Si la pendiente es mayor a 1, se fija su valor a 1*/
					if (pendiente_ab > 1){
						pendiente_ab = 1;
					}
					/* Si la pendiente es menor a -1, se fija su valor a -1*/
					else if (pendiente_ab < -1) {
						pendiente_ab = -1;
					}

					/* El angulo entre asteroides sera la ArcTan(pendiente) */
					angulo_ab = atan(pendiente_ab);

					/* Calculo de la fuerza resultante */
					fuerza = ((gravedad * asteroides[j].masa * asteroides[k].masa) / (pow(distancia_ab, 2.0)));
					/* Si la fuerza resultante es mayor a 100, se trunca a este valor */
					if (fuerza > 100) fuerza = 100;

					/* La fuerza resultante por el coseno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje x */
					asteroides[j].sum_fx += fuerza * cos(angulo_ab);
					asteroides[k].sum_fx -= fuerza * cos(angulo_ab);

					/* La fuerza resultante por el seno del angulo se agrega positivamente al asteroide j y se resta al asteroide k en el eje y */
					asteroides[j].sum_fy += fuerza * sin(angulo_ab);
					asteroides[k].sum_fy -= fuerza * sin(angulo_ab);
				}
			}
		}

		for (unsigned int q = 0; q < planetas.size(); q++) {
			/* Para el mismo asteroide "z", despues de evaluar la interaccion con el resto de asteroides, se hace lo propio
			para con todos los planetas */
			for (unsigned int z = 0; z < asteroides.size(); z++) {

				/* Calcula la distancia entre el asteroide "z" y el planeta "q" */
				distancia_ab = sqrt(pow(planetas[q].posx - asteroides[z].posx, 2.0) +
									pow(planetas[q].posy - asteroides[z].posy, 2.0));

					/* Calculo de la pendiente entre elementos */
					pendiente_ab = (planetas[q].posy - asteroides[z].posy) / (planetas[q].posx - asteroides[z].posx);

					/* Si la pendiente es mayor a 1, se fija su valor a 1*/
					if (pendiente_ab > 1){
						pendiente_ab = 1;
					}
					/* Si la pendiente es menor a -1, se fija su valor a -1*/
					else if (pendiente_ab < -1) {
						pendiente_ab = -1;
					}

					/* EL angulo entre asteroide y planeta sera el ArcTan(pendiente) */
					angulo_ab = atan(pendiente_ab);

					/* Calculo de la fuerza en el resultante entre ambos elementos */
					fuerza = ((gravedad * planetas[q].masa * asteroides[z].masa) / (pow(distancia_ab, 2.0)));
					/* Si la fuerza resultante es mayor a 100 se trunca a este valor */
					if (fuerza > 100) fuerza = 100;

					/* El producto de la fuerza por el coseno del angulo se agrega positvamente al sumatorio del asteroide "z" en el eje x */
					asteroides[z].sum_fx += fuerza * cos(angulo_ab);

					/* El producto de la fuerza por el seno del angulo se agrega positvamente al sumatorio del asteroide "z" en el eje y */
					asteroides[z].sum_fy += fuerza * sin(angulo_ab);
			}
		}

		for (unsigned int j = 0; j < asteroides.size(); j++) {

			/* Calculo de las nuevas aceleraciones */
			aceleracionx = asteroides[j].sum_fx / asteroides[j].masa;
			aceleraciony = asteroides[j].sum_fy / asteroides[j].masa;

			/* Calculo de las nuevas velocidades */
			asteroides[j].nueva_vx = asteroides[j].vx + (aceleracionx * intervalo_tiempo);
			asteroides[j].nueva_vy = asteroides[j].vy + (aceleraciony * intervalo_tiempo);

			/* Calculo de las nuevas posiciones */
			asteroides[j].nueva_posx = asteroides[j].posx + (asteroides[j].nueva_vx * intervalo_tiempo);
			asteroides[j].nueva_posy = asteroides[j].posy + (asteroides[j].nueva_vy * intervalo_tiempo);

			/* Si hubo rebote contra los bordes, se reposiciona el asteroide y modifica su velocidad */
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

		/* Para cada asteroide se actualizan las nuevas posiciones y velocidades. Tambien, todos los sumatorios se actualizan
		a 0 */
		for (unsigned int i = 0; i< asteroides.size(); i++) {

			asteroides[i].posx = asteroides[i].nueva_posx;
			asteroides[i].posy = asteroides[i].nueva_posy;
			asteroides[i].vx = asteroides[i].nueva_vx;
			asteroides[i].vy = asteroides[i].nueva_vy;
			asteroides[i].sum_fx = 0;
			asteroides[i].sum_fy = 0;
			cout << "Velocidad del asteroide [" << i << "]: " << asteroides[i].nueva_vx << ", " << asteroides[i].nueva_vy << "\n";
		}
	}

	/* Crea el fichero de salida out.txt y fija 3 decimales ale escribir en el */
	ofstream fs2("out.txt");
	fs2 << setprecision(3) << fixed;

	/* Escribe en el fichero de salida los valores finales de cada asteroide */
	for (unsigned int i = 0; i < asteroides.size(); i++) {

		fs2 << asteroides[i].posx << " " << asteroides[i].posy << " " << asteroides[i].vx << " " << asteroides[i].vy
				<< " " << asteroides[i].masa << endl;
	}

	/* Cierra el fichero out.txt */
	fs2.close();

	return 0;
}
