#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#define NUM_MASSAS 12 // Número de massas
typedef struct _molaS molaSimples;
typedef struct _massaS massaSimples;
typedef struct _mola mola;
typedef struct _massa massa;

int viewport_x = 30, viewport_y = 30;

struct _molaS {
	float k;
	float a;
};

struct _massaS {
	float m; // Massa
	double x; // Posição
	double v; // Velocidade
	double a; // Aceleração
};

const float L = 30; // Largura do "array" de massas-molas
massaSimples massas[NUM_MASSAS];
molaSimples molas[NUM_MASSAS];
float tempo = 400;
FILE *arquivo;
char *nome_arq = "saida.txt";
unsigned long num_passos;
unsigned long i;
double tol_min, tol_min_3, dt, dt2, dth, dt2h, fator_queda, f_max;
const float coef_arrasto = 0.0;
char CONTINUE;

void init();
void singleStep();
void evolve();
void desenha();
void desenhaMassa(float x, float y, float m, float lMassa);
void desenhaMola(float l, float r, float a, float k, float y, float lMassa);
void keyboard(unsigned char key, int x, int y);
void openAndInitFile();
void closeAndPlotFile();
void finalize();
void finalize2();

int main(int argc, char **argv) {

	init();
	atexit(finalize);
	//evolve();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Sistema massa-mola");

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, viewport_x, 0.0, viewport_y);

	glutDisplayFunc(desenha);
	glutKeyboardFunc(keyboard);
	glutCloseFunc(finalize2);
	glutIdleFunc(singleStep);
	CONTINUE = 1;

	glutMainLoop();
	//*/


	return 0;
}

void openAndInitFile() {
	register int i;
	arquivo = fopen(nome_arq, "w");
	fprintf(arquivo, "t\tU\tK\tE\t");
	for (i = 0; i < NUM_MASSAS; i++) {
		fprintf(arquivo, "m_%02d\t", i);
	}
	fprintf(arquivo, "m_00");
	fprintf(arquivo, "\n");
}

// Inicializa as posições, velocidades, acelerações das massas e
// disposições das molas
void init() {
	num_passos = 200 * tempo;
	tol_min = L / (NUM_MASSAS * 50);
	tol_min_3 = 3 * tol_min;
	dt = tempo / num_passos;
	dt2 = dt * dt;
	dth = dt / 2;
	dt2h = dt2 / 2;
	fator_queda = 5 / tol_min;
	f_max = 7.0;

	const float m = 1.0, k = 0.1, dx = L / NUM_MASSAS, a = 5.0;
	float p_total = 0.0, m_total = 0.0, v_cm;
	register int i = 0;

	for (i = 0; i < NUM_MASSAS; i++) {
		massas[i].m = m;
		massas[i].x = i * dx;
		massas[i].v = 0.0;
		massas[i].a = 0.0;
		molas[i].k = k;
		molas[i].a = a;

	}

	massas[0].x = 0.5 * dx;
	massas[0].v = 0.25 * dx;
	if (NUM_MASSAS > 2) {
		molas[2].k = 1.7 * k;
	}
	if (NUM_MASSAS > 7) {
		molas[7].k = 1.2 * k;
	}

	for (i = 0; i < NUM_MASSAS; i++) {
		p_total += massas[i].m * massas[i].v;
		m_total += massas[i].m;
	}

	v_cm = p_total / m_total;
	printf("%f\n", v_cm);

	p_total = 0;

	for (i = 0; i < NUM_MASSAS; i++) {
		massas[i].v -= v_cm;
		p_total += massas[i].m * massas[i].v;
	}
	v_cm = p_total / m_total;
	printf("%f\n", v_cm);


	openAndInitFile();

}

void singleStep() {
	register int j;
	float dF;
	double forcas[NUM_MASSAS];

	float U = 0, K = 0, E = 0;

	for (j = 0; j < NUM_MASSAS; j++) {
		forcas[j] = 0;
	}
	for (j = 0; j < NUM_MASSAS; j++) {
		massas[j].x += dt * massas[j].v + dt2h * massas[j].a;
	}
	for (j = 0; j < NUM_MASSAS; j++) {
		int next;
		float x = massas[j].x;
		float x_next;
		float dx;
		if (j == NUM_MASSAS - 1) {
			next = 0;
			x_next = L + massas[next].x;
		} else {
			next = j + 1;
			x_next = massas[next].x;
		}

		dx = x_next - x;

		dF = molas[j].k * (dx - molas[j].a);

		forcas[j] += dF / massas[j].m;
		forcas[next] -= dF / massas[next].m;

		// efeito de arrasto (dissipação de enrgia pelo ar)
		forcas[j] -= coef_arrasto * massas[j].v;

		U += dF * dF / (2 * molas[j].k);

		if (dx < tol_min_3) {
			float expon = exp(-fator_queda * fabs(dx)) * f_max;
			forcas[j] -= fator_queda * expon;
			forcas[next] += fator_queda * expon;

			U += expon;
		}

	}
	// Atualização das massas
	for (j = 0; j < NUM_MASSAS; j++) {
		//massas[j].x += dt*massas[j].v + dt2h*massas[j].a;
		massas[j].v += dth * (forcas[j] + massas[j].a);
		massas[j].a = forcas[j];

		K += (massas[j].m) * (massas[j].v) * (massas[j].v) * (0.5);
	}
	E = U + K;
	fprintf(arquivo, "%f\t", i * dt);
	fprintf(arquivo, "%f\t", U);
	fprintf(arquivo, "%f\t", K);
	fprintf(arquivo, "%f\t", E);
	for (j = 0; j < NUM_MASSAS; j++) {
		fprintf(arquivo, "%f\t", massas[j].x);
	}
	fprintf(arquivo, "%f", massas[0].x + L);
	fprintf(arquivo, "\n");

	if (!(i % 3)) {
		desenha();
	}
	if (CONTINUE) {
		i++;
	} else {
		closeAndPlotFile();
	}

}

void closeAndPlotFile() {
	register int i;
	char *nome_arq_plot1 = "E_molas.gnu", *nome_arq_plot2 = "x_molas.gnu";
	char comando[150] = "";

	fclose(arquivo);
	arquivo = fopen(nome_arq_plot1, "w");
	fprintf(arquivo, "set encoding iso_8859_1\n");
	fprintf(arquivo, "plot 'saida.txt' using 1:2 title 'U' with lines, \\\n");
	fprintf(arquivo, "     'saida.txt' using 1:3 title 'K' with lines, \\\n");
	fprintf(arquivo, "     'saida.txt' using 1:4 title 'E' with lines\n");
	fprintf(arquivo, "#set terminal png\n");
	fprintf(arquivo, "#set output 'energias.png'\n");
	fprintf(arquivo, "#replot\n");
	fclose(arquivo);
	strcat(comando, "gnome-terminal -e \"");
	strcat(comando, "gnuplot ");
	strcat(comando, nome_arq_plot1);
	strcat(comando, " - \" & ");
	system(comando);
	arquivo = fopen(nome_arq_plot2, "w");
	fprintf(arquivo, "plot \\\n");
	for (i = 0; i < NUM_MASSAS; i++) {
		fprintf(arquivo,
				"\t'saida.txt' using 1:%i title 'm%02d' with lines, \\\n", i
						+ 5, i);
	}
	fprintf(arquivo,
			"\t'saida.txt' using 1:%i title 'm00' with lines lt 2 lw 0.5\n", i
					+ 5);
	fprintf(arquivo, "#set terminal png\n");
	fprintf(arquivo, "#set output 'posicoes.png'\n");
	fprintf(arquivo, "#replot\n");
	fclose(arquivo);
	comando[0] = '\0';
	strcat(comando, "gnome-terminal -e \"");
	strcat(comando, "gnuplot ");
	strcat(comando, nome_arq_plot2);
	strcat(comando, " - \" & ");
	system(comando);
}

void evolve() {
	for (i = 0; i < num_passos; i++) {
		singleStep();
	}
}
void finalize() {
	printf("exit(0);");
	closeAndPlotFile();
	exit(0);
}

void finalize2() {
	printf("exit(2);");
	CONTINUE = 0;

}

void desenha(void) {
	int i;
	float y;
	float lMassa = L / (NUM_MASSAS * 15); // Proporcional à distância de interação do "potencial de contato"
	//float k;
	//float x_c, y_c, r, l, ang, a;
	glClear(GL_COLOR_BUFFER_BIT);

	y = 12;
	glColor4f(0.5, 0.5, 1.0, 0.5);

	glPushMatrix();
	glScalef(viewport_x / L, 1.0, 1.0);
	glColor4f(0.3, 0.3, 0.3, 0.5);
	desenhaMassa(massas[0].x, y, massas[0].m, lMassa);
	glColor4f(0.5, 0.5, 1.0, 0.5);
	for (i = 1; i < NUM_MASSAS; i++) {
		desenhaMassa(massas[i].x, y, massas[i].m, lMassa);
	}
	glColor4f(0.3, 0.3, 0.3, 0.5);
	desenhaMassa(massas[0].x + L, y, massas[0].m, lMassa);

	glColor3f(1.0, 0.5, 0.5);
	for (i = 0; i < NUM_MASSAS - 1; i++) {
		desenhaMola(massas[i].x, massas[i + 1].x, molas[i].a, molas[i].k, y,
				lMassa);
	}
	desenhaMola(massas[i].x, massas[0].x + L, molas[i].a, molas[i].k, y, lMassa);
	desenhaMola(massas[i].x - L, massas[0].x, molas[i].a, molas[i].k, y, lMassa);
	if (NUM_MASSAS > 1) {
		desenhaMola(massas[0].x + L, massas[1].x + L, molas[0].a, molas[0].k,
				y, lMassa);
	}

	glPopMatrix();
	/*
	 glColor4f(0.0, 0.0, 0.0, 0.0);
	 glBegin(GL_POINTS);
	 for (i = 1; i <= viewport_x; i++) // Desenha os pontos das coordenadas
	 int j;
	 for (j = 1; j <= viewport_y; j++)
	 glVertex2i(i, j);
	 glEnd();
	 */
	glFlush();

}

void desenhaMassa(float x, float y, float m, float lMassa) {
	glRectf(x - lMassa, y, x + lMassa, y + m);
}

void desenhaMola(float l, float r, float a, float k, float y, float lMassa) {
	float dxMola = 0.5 * lMassa;
	float w, ang, rad;
	float xr = l + lMassa + dxMola, xl = r - lMassa - dxMola;
	float Dx = xl - xr;

	//glRectf(l + lMassa, y, r - lMassa, y + 0.5);


	ang = (2 * a) * M_PI;
	rad = 0.5;

	glBegin(GL_LINE_STRIP);

	glVertex2f(l + lMassa, y + 0.5);
	glVertex2f(l + lMassa + dxMola, y + 0.5);

	for (w = 0; w <= ang; w += 0.1) {
		float x = xr + Dx * w / ang;
		glVertex2f(x, y + 0.5 + rad * sin(w));
	}

	glVertex2f(r - lMassa - dxMola, y + 0.5);
	glVertex2f(r - lMassa, y + 0.5);

	glEnd();

}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		finalize();
	};
	glutPostRedisplay();
}

struct _massa {
	float m;
	float x;
	mola* mola_l;
	mola* mola_r;
};

struct _mola {
	float k;
	float a;
	massa* massa_l;
	massa* massa_r;
};
