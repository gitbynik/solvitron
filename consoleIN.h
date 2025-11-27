#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#ifndef _INH_
#define _INH_
#include "dsGraph.h"
#include "tree.h"
#include "createMat.h"

extern int maxargs;
extern bool running;

void plsseparate(char*, char(*)[10]);
void createnode(int);
bool createComp(unsigned short int, int, int, double);
void handleinput(char*, char(*)[10],gsl_matrix*);
void test(graphNode);
void plot_row(gsl_matrix *,double,int,int); 
enum components : unsigned short int{Resistor = 1, Capacitor = 1 <<1, Inductor = 1 <<2, VoltSource = 1 <<3};
#endif