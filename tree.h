#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifndef _TREEH_
#define _TREEH_
#include "dsGraph.h"

extern int graphNodecount ;
extern graphNode graphHead;
extern graphNode graphTail;
extern graphNode cotreeHead;
extern graphNode treeRoot;

void DFStree(graphNode);
void push(graphNode);
graphNode pop();
graphNode ltog(linkedNode);

#endif