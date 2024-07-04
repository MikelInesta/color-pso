/*
	Miguel Iñesta Garza
	Sistemas Inteligentes
	Compresión de color en imágenes con PSO
	2022
*/

/*---------------------------------Librerías-------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*-----------Constantes y variables globales para pso colores-------------------------*/


	//Numero de pixeles de la imagen
	#define N_PUNTOS 262144
	 
	//Numero de columnas de pixeles de la imagen
	int n_c;
	
	//Numero de filas de pixeles de la imagen
	int n_f;
	
	/*Numero de colores de la paleta cuantizada 8 bits de profundidad*/
	#define R 256
	
	/*Tamaño de la poblacion de particulas*/
	#define P 15
		
	/*Valores para posicion random*/
	#define lim_min 0 		
	#define lim_max 255
	
	/*Valores para velocidad random*/
	#define v_min -5
	#define v_max 5
		
	/*Inercia*/
    #define paramW  0.729 

	/* Peso componente cognitiva*/
	#define f1  2.05 
        
	/* Peso componente social */
	#define f2  2.05 
		
	/* numero maximo de iteraciones de un algoritmo enjambre */
	int max_iterations = 20;
	
	//Vectores para almacenar colores de cada pixel de la imagen original
	int ro[N_PUNTOS], ve[N_PUNTOS], az[N_PUNTOS];

    int  ro_cuantizada[N_PUNTOS], ve_cuantizada[N_PUNTOS], az_cuantizada[N_PUNTOS];

/*------------------------------REGISTROS--------------------------------------------*/

	//Registro que contiene la informacion de una particula
	struct particle {

	//posicion actual (color (1-8) y componentes de cada color(0-255))
	int x[R][3];

	//velocidad 
	double v[R][3];

	int b[R][3];

};

typedef struct particle PARTICLE;

/*------------------------------PROTOTIPOS-------------------------------------------*/

double funcion_mse();
void escribir_pixels_imagen(char nombre_fich[]);
double distancia_euclidea(int p1r, int p2r, int p1g, int p2g, int p1b, int p2b);
void generar_imagen_cuantizada(int paleta[R][3]);
void leer_pixels_imagen(char nombre_fich[]);
void spawn_swarm(PARTICLE particle_array[], double fitness[P]);
double get_random_double(double min, double max);
void fitnessCheck(PARTICLE particle_array[], double fitness[]);
void updatePersonalBest(PARTICLE particle_array[], double fitness[]);
void updateGlobalBest(PARTICLE particle_array[], int g[R][3], double* fitness_g);
double globalBestInit(PARTICLE particle_array[], double fitness[], int g[R][3]);
void update_speed_position(PARTICLE particle_array[], int g[R][3]);
void show_particle_fitness(double fitness[]);
void show_PSO_result(int g[R][3], double fitness_g);
void PSO();

/*----------------------------------MAIN----------------------------------------------*/

void main() {
    
    srand48(time(NULL));
    srand(time(NULL));
    PSO();
    
}

/*------------------------------FUNCION PSO-------------------------------------------*/
void PSO() {

	//leemos las componentes de cada pixel de una imagen y los almacenamos en las variables globales ro[] ve[] az[]
    leer_pixels_imagen("datos_mandril.txt");

    // vector para almacenar las P particulas del enjambre
    PARTICLE particle_array[P];

    // matriz con la mejor solucion encontrada por el enjambre hasta ahora (es una paleta de 8 colores y 24 componentes)
    int g[R][3];

    // fitness de cada particula del enjambre
    double fitness[P];

    //fitness de la mejor solucion global que esta almacenada en el vector g (real que devuelve mse)
    double fitness_g;

    // contador de iteraciones del algoritmo
    int i;

    //1.- inicializar la poblacion de partículas
    spawn_swarm(particle_array, fitness);
	
	//AQUI ESTA EL PRIMER ERROR ------------------------------------------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!
    //g[] no tiene valor inicial asi que añadimos una funcion para definirlo
    fitness_g = globalBestInit(particle_array, fitness, g);


    //2.- bucle del pso, cada iteracion es un paso de mejora
    for (i = 0; i < max_iterations; i++) {


        /*En la primera iteracion de PSO no aplico esta operacion, pues acabo de 
        calcular los fitness al determinar el valor inicial del vector g */
        if (i != 0) {

            //2-a evaluar el fitness de cada particula
            fitnessCheck(particle_array, fitness);

            //2-b actualizar la mejor solucion personal de cada particula
            updatePersonalBest(particle_array, fitness);

            //2-c actualizar la mejor solucion global
            updateGlobalBest(particle_array, g, &fitness_g); 

        }

        //2-d actualizar la velocidad y la posicion de cada particula
        update_speed_position(particle_array, g);

        //Muestro la mejor solucion hasta el momento
        printf("\nIteracion %d ---------> ", i+1);
        show_PSO_result(g, fitness_g);

    }

    //Muestro la mejor solucion final
    printf("\n<------------- MEJOR SOLUCION ENCONTRADA POR PSO ------------->");
    show_PSO_result(g, fitness_g);
    printf("\n");

    escribir_pixels_imagen("mandril_cuantizado.txt");

}

/*---------------------FUNCION MOSTRAR RESULTADO---------------------------------------//
		funcion que muestra el resultado obteido al aplicar el algoritmo PSO
		params:
		    g[][]: mejor resultado encontrado por el enjambre (matriz)
		    fitness_g : la calidad del mejor resultado encontrado
		return:
		    void
		[[ probablemente no tenga errores ]]
*/
void show_PSO_result(int g[R][3], double fitness_g) {
     
    /*
    printf("\n(");
    for(int i=0; i<R; i++){
		printf("\n(");
	    printf("%d ", g[i][0]);
	    printf("%d ", g[i][1]);
	    printf("%d ", g[i][2]);
	    printf(")");
    }
    printf(")");
    */
    printf(" FITNESS-> %lf", fitness_g);
}

/*-----------------FUNCION MOSTRAR FITNESS DE TODAS LAS PARTICULAS------------------*/
// [[ probablemente no tenga errores  ]]
void show_particle_fitness(double fitness[]) {
    int i;
    for (i = 0; i < P; i++) {
        printf("Par->%lf ", fitness[i]); 
    }
}

/*----------FUNCION QUE ACTUALIZA VELOCIDAD Y POSICION DE CADA PARTICULA--------*/
// [[ PARECE NO TENER ERRORES ]]
void update_speed_position(PARTICLE particle_array[], int g[R][3]) {
    int i, j, k;

    for (i = 0; i < P; i++) {
        for (j = 0; j < R; j++) {
            for(k=0; k<3; k++){
                particle_array[i].v[j][k] = paramW * particle_array[i].v[j][k] + f1 * drand48() * 
                                            (particle_array[i].b[j][k] - particle_array[i].x[j][k]) +
                                            f2 * drand48() * (g[j][k] - particle_array[i].x[j][k]);
                      
                //Comprobamos que la velocidad no se salga del rango                      
		        if(particle_array[i].v[j][k]<v_min){
		        	particle_array[i].v[j][k]=v_min;
		        }  
		        if(particle_array[i].v[j][k]>v_max){
		        	particle_array[i].v[j][k]=v_min;
		        }

                particle_array[i].x[j][k] += particle_array[i].v[j][k];

                //comprobamos que la posicion no se salga del rango
	  	     	if(particle_array[i].x[j][k]<lim_min){
			       	particle_array[i].x[j][k]=lim_min;
		        }  
		        if(particle_array[i].x[j][k]>lim_max){
			       	particle_array[i].x[j][k]=lim_min;
	            }

            }
        }
    }
}


/*---------------------FUNCION QUE INICIALIZA G--------------------------------*/
double globalBestInit(PARTICLE particle_array[], double fitness[], int g[R][3]) {

	//Contadores
    int i,j;
    
    //Mejor fitness encontrado hasta ahora
    double bestFit;
    
    //La posicion de la particula que me ha dado el mejor fitness hasta ahora
    int bestParticle;

    //tomo una particula como mejor inicial (la primera)
    //posicion en el vector de particulas
    bestParticle = 0;
    
    

    //El fitness del vector b de la primera particula
    generar_imagen_cuantizada(particle_array[0].b);
    bestFit = fitness[0] = funcion_mse();

    //Comparo todas las particulas para ver si hay una mejor que la inicial
    //Empiezo en 1 porque la 0 la tomé como inicial
    for (i = 1; i < P; i++) {
        //calculo el fitness de cada particula
        generar_imagen_cuantizada(particle_array[i].b);
        fitness[i] = funcion_mse();

        //Si es mejor que todos los fitness calculados la asignamos como 
        if (fitness[i] < bestFit) {
            //actualizo la mejor particula
            bestParticle = i;
            //copio su fitness
            bestFit = fitness[i];
        }
    }
    //Se que la particula con mejor posicion personal es la identificada por 'bestParticle'
    //copio en g el vector b de dicha particula
    for (i = 0; i < R; i++) {
        for(j=0; j< 3; j++){
            g[i][j] = particle_array[bestParticle].b[i][j];
        }
    }
    return bestFit;

}

void updateGlobalBest(PARTICLE particle_array[], int g[R][3], double* fitness_g) { 

    int i, j; //contador de particulas

    double a_fitness; //fitness de una posicion (de la mejor personal de alguna particula)

    int bestParticle = -1; //particula cuyo vector b es mejor que el vector g actual
    //Le doy un valor inicial para indicar que no se encontro mejor

    for (i = 0; i < P; i++) {

        //calculo el fitness de su mejor posicion personal
        generar_imagen_cuantizada(particle_array[i].b);
        a_fitness = funcion_mse();

        //Si es mejor que el fitness de g
        if (a_fitness < *fitness_g) {

            //anoto el indice de la particula cuyo b tiene el mejor fitness hasta ahora
            bestParticle = i;

            //copio su fitness
            *fitness_g = a_fitness;

            //Esto lo tengo que mejorar!!!!!!! Pero funciona

        }
    }
    if (bestParticle != -1) {
        //actualizo el valor a g (la mejor posicion personal de la partícula i)
        for (i = 0; i < R; i++) {
            for(j=0; j< 3; j++){
                g[i][j] = particle_array[bestParticle].b[i][j];
            }
        }
    }
}

//Funcion que actualiza la mejor posicion global de cada particula
//comparo el fitness actual de (x) con el mejor fitness de la mejor posicion personal (b).
//Si el fitness de x es mejor lo asigno a b
void updatePersonalBest(PARTICLE particle_array[], double fitness[]) {
    for (int i = 0; i < P; i++) {
        //comparo el fitness de x y b
        //Si el fitness de x es menor, copio x sobre b, esto se debe a que el problema busca minimizacion
        generar_imagen_cuantizada(particle_array[i].b);
        if (fitness[i] < funcion_mse()) {
            for (int j = 0; j < R; j++) {
                for(int k=0; k<3; k++){
                    particle_array[i].b[j][k] = particle_array[i].x[j][k];
                }
            }
        }
    }
}


/*Calcular el fitness de cada particula*/
void fitnessCheck(PARTICLE particle_array[], double fitness[]) {
    
    //Determinar que color de la paleta es mas parecido a cada pixel de la imagen original
    for(int i=0; i<P; i++){
            generar_imagen_cuantizada(particle_array[i].x);
            fitness[i] = funcion_mse();
        }
}

double get_random_double(double min, double max) {
    //srand(time(NULL)); // MIM: lo quité porque lo estaba inicializando igual, con lo que daba el mismo rand() :(
    float scale = rand() / (float)RAND_MAX;
    return min + scale * (max - min); // MIM: cambie el órden, porque siempre da negativo así (max es más grande que min).
}

//inicializa el enjambre de particulas
void spawn_swarm(PARTICLE particle_array[], double fitness[]) {

    for (int i = 0; i < P; i++) {
        for (int j = 0; j < R; j++) { 
        	for(int k=0; k< 3; k++){
        		particle_array[i].x[j][k] = get_random_double(lim_min, lim_max);
        		particle_array[i].v[j][k] = get_random_double(v_min, v_max);
        		particle_array[i].b[j][k] = particle_array[i].x[j][k];
        	}
        }
        generar_imagen_cuantizada(particle_array[i].x);
        fitness[i]=funcion_mse();
    }
}

/*-------------------------------COSAS DE PSO--------------------------------------*/



/*----------------------------------------------------------------------------------/
/-----------------------------------COSAS DE IMAGEN---------------------------------/
/----------------------------------------------------------------------------------*/


/*----------------------FUNCION PARA LEER PIXELES DE UNA IMAGEN--------------------*/
void leer_pixels_imagen(char nombre_fich[]){
	
	FILE *fd; /* file handler*/
	int i;	/*contador*/

	if((fd = fopen(nombre_fich, "r"))==NULL){
		printf("Error al abrir el fichero \n");
		exit(1);
	}

	fscanf(fd, "%d %d", &n_f, &n_c); //la primera linea indica el numero de filas y columnas de la matriz

	int N = n_f * n_c;

	if (N > N_PUNTOS){

		printf("\n ERROR AL LEER: hay demasiados puntos. Incrementar la ct. N_PUNTOS y recompilar");
		fclose(fd);
		exit(1);
	}
    else{
        for(i=0; i<N; i++){
            fscanf(fd, "%d %d %d", &ro[i], &ve[i], &az[i]);

        }
        fclose(fd);
    }

}

/*
	Funcion que genera la imagen cuantizada
*/
void generar_imagen_cuantizada(int paleta[R][3]){

	int i, j; //contador
	int mejor_color_paleta; //indice del color mas parecido de la paleta (entre 0 y r-1)

	double min_distancia, //distancia entre un pixel de la imagen y el color mas parecido de la paleta
			dist; //una distancia

    

	//para cada punto de la imagen
	for(i=0; i<N_PUNTOS; i++){
		//comparo el pixel i-esimo de la imagen original con todos los colores de la paleta y me quedo con el mas parecido

        mejor_color_paleta = 0;
		min_distancia = distancia_euclidea(ro[i] ,paleta[0][0], ve[i], paleta[0][1], az[i], paleta[0][2]);

		// miro los r-1 colores restantes de la paleta
		for(j=1; j<R; j++){

			//calculo la distancia entre el pixel i-esimo y el color j-esimo
			dist = distancia_euclidea(ro[i], paleta[j][0], ve[i], paleta[j][1], az[i], paleta[j][2]);

			//comparo con la mejor distancia que tengo almacenada
			if (dist < min_distancia){
				/*he encontrado un color de la paleta mas parecido al pixel i-esimo
				de la imagen original asi que lo tomo como nuevo mejor color
				*/
				mejor_color_paleta = j;
				min_distancia = dist;
			}

		}
		// en este punto conozco el color de la paleta mas parecido al pixel i-esimo de la imagen original
		
		// ahora defino el color del pixel que ocupara la posicion i en la imagen cuantizada
		ro_cuantizada[i] = paleta[mejor_color_paleta][0];
		ve_cuantizada[i] = paleta[mejor_color_paleta][1];
		az_cuantizada[i] = paleta[mejor_color_paleta][2];
        
	}
}

//--------------------Distancia euclidea entre 3 pares de enteros-----------------------------
double distancia_euclidea(int p1R, int p2R, int p1G, int p2G, int p1B, int p2B){

	return sqrt( pow(p1R-p2R,2) + pow(p1G - p2G,2) + pow(p1B-p2B,2));

}

void escribir_pixels_imagen(char nombre_fich[]){

	FILE *fd; // descriptor del fichero
	int i; //contador

	if((fd = fopen(nombre_fich, "w"))==NULL){

		printf("Error al abrir el fichero\n");
		exit(-2);
	}

	fprintf(fd, "%d %d\n", n_f, n_c); //la primera linea indica el numero de filas y columnas de la imagen que voy a volcar al disco

	//Se escribe el color de los sucesivos pixels
	// Escribo 1 pixel en cada linea
	for(i=0; i<N_PUNTOS; i++){

		fprintf(fd, "%d %d %d\n", ro_cuantizada[i], ve_cuantizada[i], az_cuantizada[i]);
		
	}
    fclose(fd);
}

//-----------------------FUNCION OBJETIVO (ERROR DE CUANTIFICACION)---------------------/
// [[ PARECE NO TENER ERRORES ]]
double funcion_mse(){
	int i;
	double mse=0;

    for(i = 0; i<N_PUNTOS; i++){
        mse+=pow((ro[i]-ro_cuantizada[i]),2)+pow((ve[i]-ve_cuantizada[i]),2)+pow((az[i]-az_cuantizada[i]),2);
    }

    mse/=N_PUNTOS;

    return mse;
	
}	

