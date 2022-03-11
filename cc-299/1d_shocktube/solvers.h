#ifndef SOLVERS_H
#define SOLVERS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"

void centeredScheme(results *, grid *);
void laxWendroff(results *, grid *);
void macCormack(results *, grid *);
void stegerWarming(results *, grid *);
void vanLeerNonMUSCL(results *, grid *);
void liouAUSMplus(results *, grid *);
void roeMethod(results *, grid *);
void hartenTVD(results *, grid *);
void exactSolution(results *, grid *);

#endif
