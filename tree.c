#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "tree.h"

graphNode stackptr = NULL;

void push(graphNode topush)
{
	if(topush->nextinstack != NULL);
	topush->nextinstack = stackptr;
	topush->visited = true;
	stackptr = topush;
}

graphNode pop()
{
	if(stackptr == NULL) return NULL;
	graphNode tmp = stackptr;
	stackptr = stackptr->nextinstack;
	tmp->nextinstack = NULL;
	return tmp;
}

graphNode ltog(linkedNode link)
{
	int x;
	if(link->isforwardlink)	x = link->destindex;
	else	x = link->hostindex;
	graphNode tmp = searchgraphNode(x);
	return tmp;
}

void DFStree(graphNode start)
{
	push(start);
	linkedNode templink = start->nextlink;
	while(templink!=NULL)
	{
		if(ltog(templink)->visited)
		{
			int x = (start->index == templink->destindex) ? templink->hostindex : templink->destindex;
			if(stackptr->nextinstack!=searchgraphNode(x))	templink->iscoTreeLink = true;
		}
		else	DFStree(ltog(templink));
		templink=templink->nextlink;
	}
	pop();
}