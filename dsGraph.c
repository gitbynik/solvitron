#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "dsGraph.h"
//deletegraphnode does not edit link nodes, plis fix
graphNode graphTail = NULL;
graphNode cotreeHead = NULL;
graphNode treeRoot = NULL;
extern graphNode graphHead;
extern int graphNodecount;
extern int edgescount;

void freeall()
{
	graphNode tmp = graphHead;
	graphNode tmp2 = tmp;
	while(tmp!=NULL)
	{
		linkedNode tmplink = tmp->nextlink;
		linkedNode tmplink2 = tmplink;
		while(tmplink!=NULL)
		{
			tmplink2 = tmplink2->nextlink;
			free(tmplink);
			tmp->nextlink=NULL;
			tmplink = tmplink2;
		}
		tmp2 = tmp2->next;
		free(tmp);
		tmp = tmp2;
	}
	graphNodecount =0;
	graphHead = NULL;
	graphTail = NULL;
}

graphNode creategraphNode(int index)
{
	graphNode extratempo = searchgraphNode(index);
	if(extratempo!=NULL) return extratempo; //prevents duplication of nodes
	graphNode tmpnode = (graphNode) malloc(sizeof(struct plsgraphNode));
	graphNodecount++;
	tmpnode->index = index;
	tmpnode->nextlink = NULL;
	tmpnode->visited = false;
	tmpnode->nextinstack = NULL;
	if(graphHead!=NULL)
	{
		graphHead->prev = tmpnode;
		tmpnode->next=graphHead;
	}
	else tmpnode->next = NULL;
	graphHead = tmpnode;
	tmpnode->prev = NULL;
	if(graphTail==NULL) graphTail = tmpnode;
	return tmpnode;
}

void deletegraphNode(graphNode tmpnode)
{
	graphNodecount--;
	if(graphHead==graphTail)
	{
		graphHead = NULL;
		graphTail = NULL;
	}
	else if(tmpnode->prev==NULL)
	{
		graphHead = tmpnode->next;
		graphHead->prev = NULL;
	}
	else if(tmpnode->next == NULL)
	{
		graphTail = tmpnode->prev;
		graphTail->next =NULL;
	}
	else
	{
		(tmpnode->prev)->next = tmpnode->next;
		(tmpnode->next)->prev = tmpnode->prev;
	}
	free(tmpnode);
}

graphNode searchgraphNode(int index) //search by index
{
	graphNode tmpnode = graphHead;
	while(tmpnode != NULL)
	{
		if(tmpnode->index==index) return tmpnode;
		tmpnode = tmpnode->next;
	}
	return NULL;
}

void addEdge(unsigned short int comp, double value, graphNode source, graphNode dest)
{
	edgescount++;
	linkedNode result = (linkedNode) malloc(sizeof(struct plslinkedNode));
	result->nextlink = NULL;
	result->prevlink = NULL;
	result->comp = comp;
	result->branchIndex = 0;
	result->value = value;
	result->destindex = dest->index;
	result->hostindex = source->index;
	result->iscoTreeLink = false;
	result->isforwardlink = true;

	if(source->nextlink==NULL) source->nextlink = result;
	else
	{
		linkedNode tmp = source->nextlink;
		while(tmp->nextlink!=NULL)
		{
			tmp=tmp->nextlink;
		}
		tmp->nextlink = result;
		result->prevlink=tmp;
	}

	linkedNode undirResult = (linkedNode) malloc(sizeof(struct plslinkedNode));
	undirResult->nextlink = NULL;
	undirResult->prevlink = NULL;
	undirResult->comp = comp;
	undirResult->branchIndex = 0;
	undirResult->value = value;
	undirResult->destindex = dest->index;//plscheck
	undirResult->hostindex = source->index;//plscheck
	undirResult->iscoTreeLink = false;
	undirResult->isforwardlink = false;
	
	if(dest->nextlink==NULL) dest->nextlink = undirResult;
	else
	{
		linkedNode tmp = dest->nextlink;
		while(tmp->nextlink!=NULL)
		{
			tmp=tmp->nextlink;
		}
		tmp->nextlink = undirResult;
		undirResult->prevlink = tmp;
	}
}

void removeEdge(graphNode source, graphNode dest)
{
	edgescount--;
	int dstindex = dest->index;
	int srcindex = source->index;
	linkedNode tmplink = source->nextlink;
	linkedNode tmplink2 = tmplink;
	if(tmplink->destindex==dstindex) source->nextlink=tmplink->nextlink;
	else
	{
		while(tmplink->destindex!=dstindex)
		{
			tmplink2 = tmplink;
			tmplink = tmplink->nextlink;
		}
		tmplink2->nextlink = tmplink->nextlink;
	}
	free(tmplink);

	tmplink = dest->nextlink;
	tmplink2 = tmplink;
	if(tmplink->hostindex==srcindex) dest->nextlink=tmplink->nextlink;
	else
	{
		while(tmplink->hostindex!=srcindex)
		{
			tmplink2 = tmplink;
			tmplink = tmplink->nextlink;
		}
		tmplink2->nextlink = tmplink->nextlink;
	}
	free(tmplink);
}

void printGraph(graphNode start)
{
	graphNode tmp = start;
	while(tmp!=NULL)
	{
		linkedNode tmplink = tmp->nextlink;
		printf("\nNode %d",tmp->index);
		while(tmplink!=NULL)
		{
			if(tmplink->isforwardlink)	printf(" is connected to %d",tmplink->destindex);
			else	printf(" is connected from %d",tmplink->hostindex);
			if(tmplink->iscoTreeLink)	printf(" (cotree) ");
			tmplink=tmplink->nextlink;
		}
		tmp=tmp->next;
	}
}

void inigraph(graphNode graphHead)
{
	graphNode temp = graphHead;
	while(temp!=NULL)
	{
		linkedNode tmplink = temp->nextlink;
		while(tmplink!=NULL)
		{
			tmplink->branchIndex = 0;
			tmplink = tmplink->nextlink;
		}
		temp=temp->next;
	}
}