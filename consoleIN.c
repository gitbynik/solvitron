#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "consoleIN.h"

graphNode graphHead = NULL;
int graphNodecount = 0;
int edgescount = 0;
int gnd=1;

void plot_row(gsl_matrix *A,double timewindow, int numofentries, int rowtoplot) 
{
	FILE *gp = popen("gnuplot -persistent", "w");

	fprintf(gp, "set title 'current in branch no. %d'\n", 1+rowtoplot);
	fprintf(gp, "set xlabel 'time'\n");
	fprintf(gp, "set ylabel 'current in mA'\n");
	fprintf(gp, "plot '-' using 1:2 with lines lw 2 title 'Row %d'\n", rowtoplot);

	for (int j = 0; j < numofentries; j++) 
	{
		fprintf(gp, "%lf %lf\n", j*timewindow/numofentries, (1000*((gsl_matrix_get(A,rowtoplot,j)))));
	}

	fprintf(gp, "e\n");
	fflush(gp);
	pclose(gp);
}

void plsseparate(char* in, char (*delimout)[10])
{
	char* token = strtok(in, " ");
	for(int i=0;i<maxargs;i++)
	{
		strncpy(delimout[i], token, 9);
		delimout[i][9] = '\0';
		token = strtok(NULL, " ");
		if(token == NULL) return;
	}
}

void createnode(int nodeindex)
{
	if(searchgraphNode(nodeindex)!=NULL)
	{
		printf("Node already exists");
		return;
	}
	else if(nodeindex == 0)	
	{
		printf("Node index not entered or incorrect, please try again\n");
		return;
	}
	creategraphNode(nodeindex);
	printf("Created Node with index %d \n", nodeindex);
}

bool createComp(unsigned short int comp, int nodeindex1, int nodeindex2, double value)
{
	graphNode node1 = searchgraphNode(nodeindex1);
	graphNode node2 = searchgraphNode(nodeindex2);
	if(comp == 69) return false;
	if(comp == 0 || value == 0 || nodeindex1 == 0 || nodeindex2 == 0 || node1 == NULL || node2 == NULL)
	{
		printf("Incorrect node index/ comp value, please try again\n");
		return false;
	}
	printf("Added component of %f value from node %d to node %d \n", value,nodeindex1,nodeindex2);
	addEdge(comp,value,node1,node2);
	return true;
}

void handleinput(char* input, char (*din)[10],gsl_matrix* ib)
{
	input[strcspn(input, "\n")] = 0;
	memset(din, 0, sizeof(din));
	plsseparate(input,din);
	unsigned short int comp = 69;
	if (strcmp(input, "exit") == 0) 
	{
		printf("Exiting program.\n");
		running = false;
	}
	else if(strcmp(input, "print") == 0) printGraph(graphHead);
	else if(strcmp(input, "tree") == 0) DFStree(graphHead);
	else if(strcmp(input, "clear") == 0) printf("\e[1;1H\e[2J");
	else if(strcmp(input, "refnode") == 0) gnd = atoi(din[1]);
	else if(strcmp(din[0],"res") == 0) comp = Resistor;
	else if(strcmp(din[0],"cap") == 0) comp = Capacitor;
	else if(strcmp(din[0],"ind") == 0) comp = Inductor;
	else if(strcmp(din[0],"volt") == 0) comp = VoltSource;
	else if(strcmp(din[0],"node") == 0) createnode(atoi(din[1]));
	else if(strcmp(din[0],"plot") == 0) plot_row(ib,0.3,1000,atoi(din[1]));
	if(!createComp(comp, atoi(din[2]), atoi(din[3]), strtod(din[1],NULL))) return;
}

void test(graphNode graphHead)
{
	createnode(1);
	createnode(2);
	createnode(3);
	createnode(4);
	gnd = 4;
	createComp(Resistor, 1, 2, 100);
	createComp(Inductor, 3, 4, 1);
	createComp(Resistor, 2, 3, 100);
	createComp(Capacitor, 2, 3, 0.00001);
	createComp(VoltSource, 4, 1, 10);
	printGraph(graphHead);
	printf("\e[1;1H\e[2J");
}