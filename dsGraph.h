#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifndef _GRAPH_
#define _GRAPH_

struct plslinkedNode
{
    int destindex; //which graphnode this represents, dest graphnode
    int hostindex; //source graphnode
    unsigned short int comp; //res or cap or etc
    double value; //value of resistance, capacitance etc
    int branchIndex;
    bool iscoTreeLink;
    bool isforwardlink; //true if link is hostnode->destnode
    struct plslinkedNode* nextlink;
    struct plslinkedNode* prevlink; //used for traversal inside links list
};
typedef struct plslinkedNode* linkedNode; /*these are the fictiotious nodes which are part of adjacency linked list
but since they are represent edges, these are actually the RLC components in circuit*/

struct plsgraphNode
{
    int index;
    bool visited;
    struct plsgraphNode* next;
    struct plsgraphNode* prev;
    struct plslinkedNode* nextlink;
    struct plsgraphNode* nextinstack;//used in tree making, not required for graphs
};
typedef struct plsgraphNode* graphNode; //these are the actual electrical nodes

extern graphNode graphHead;

graphNode creategraphNode(int);
void deletegraphNode(graphNode);
graphNode searchgraphNode(int);
void addEdge(unsigned short int, double,graphNode,graphNode); /*first is from which, so we go to its nextlink , second is to which, 
so destindex = this ones index*/
void removeEdge(graphNode,graphNode);
void printGraph(graphNode);
void inigraph(graphNode);
//graphNode checkConnected(graphNode,graphNode); //if first->second
void freeall();

#endif