#ifndef SOLVERS_H
#define SOLVERS_H

void centeredScheme(results *, grid *);
void laxWendroff(results *, grid *);
void macCormack(results *, grid *);
void stegerWarming(results *, grid *);
void vanLeerNonMUSCL(results *, grid *);
void liouAUSMplus(results *, grid *);
void roeMethod(results *, grid *);
void exactSolution(results *, grid *);

#endif
