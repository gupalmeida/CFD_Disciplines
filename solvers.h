#ifndef SOLVERS_H
#define SOLVERS_H

void centeredScheme(results *, grid *);
void laxWendroff(results *, grid *);
void macCormack(results *, grid *);
void stegerWarming(results *, grid *);
void exactSolution(results *, grid *);

#endif
