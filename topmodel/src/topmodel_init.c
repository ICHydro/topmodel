#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void c_flowlength(void *, void *, void *, void *, void *, void *);
extern void c_infiltration(void *, void *, void *, void *);
extern void c_sinkfill(void *, void *, void *, void *, void *, void *);
extern void c_streamorder(void *, void *, void *, void *, void *, void *);
extern void c_subcatch(void *, void *, void *, void *, void *, void *);
extern void c_topidx(void *, void *, void *, void *, void *, void *, void *);
extern void c_topmodel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void findrivers(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"c_flowlength",   (DL_FUNC) &c_flowlength,    6},
    {"c_infiltration", (DL_FUNC) &c_infiltration,  4},
    {"c_sinkfill",     (DL_FUNC) &c_sinkfill,      6},
    {"c_streamorder",  (DL_FUNC) &c_streamorder,   6},
    {"c_subcatch",     (DL_FUNC) &c_subcatch,      6},
    {"c_topidx",       (DL_FUNC) &c_topidx,        7},
    {"c_topmodel",     (DL_FUNC) &c_topmodel,     12},
    {"findrivers",     (DL_FUNC) &findrivers,      9},
    {NULL, NULL, 0}
};

void R_init_topmodel(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
