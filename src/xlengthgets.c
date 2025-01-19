/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2023  The R Core Team
 *  Copyright (C) 1995-1998  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *  
 *  Copyright (C) 2025  Merlise Clyde
 *  modified by Merlise Clyde to allow for initial values when growing vectors and
 *  better memory allocation for large vectors
 
 */


/* 
*  from R-source/builtin.c  based on xlengthgets() 
*  SEXP xlengthgets(SEXP x, R_xlen_t len)
*
*/
#include "bas.h"

SEXP xgrowvector(SEXP x, R_xlen_t len)
//  SEXP xlengthgets(SEXP x, R_xlen_t len)
  {
    R_xlen_t lenx, i;
    SEXP rval, names, xnames, t;
    if (!isVector(x) && !isList(x))
      error(_("cannot set length of non-(vector or list)"));
    if (len < 0) error(_("invalid value")); // e.g. -999 from asVecSize()
    if (isNull(x) && len > 0)
      warning(_("length of NULL cannot be changed"));
    lenx = xlength(x);
    if (lenx == len)
      return (x);
    PROTECT(rval = allocVector(TYPEOF(x), len));
    PROTECT(xnames = getAttrib(x, R_NamesSymbol));
    if (xnames != R_NilValue)
      names = allocVector(STRSXP, len);
    else names = R_NilValue;	/*- just for -Wall --- should we do this ? */
switch (TYPEOF(x)) {
case NILSXP:
  break;
case LGLSXP:
case INTSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      INTEGER(rval)[i] = INTEGER(x)[i];
    }
    else
      INTEGER(rval)[i] = 0;
    break;
case REALSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      REAL(rval)[i] = REAL(x)[i];
      if (xnames != R_NilValue)
        SET_STRING_ELT(names, i, STRING_ELT(xnames, i));
    }
    else
      REAL(rval)[i] = NA_REAL;
    break;
case CPLXSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      COMPLEX(rval)[i] = COMPLEX(x)[i];
      if (xnames != R_NilValue)
        SET_STRING_ELT(names, i, STRING_ELT(xnames, i));
    }
    else {
      COMPLEX(rval)[i].r = NA_REAL;
      COMPLEX(rval)[i].i = NA_REAL;
    }
    break;
case STRSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      SET_STRING_ELT(rval, i, STRING_ELT(x, i));
      if (xnames != R_NilValue)
        SET_STRING_ELT(names, i, STRING_ELT(xnames, i));
    }
    else
      SET_STRING_ELT(rval, i, NA_STRING);
    break;
case LISTSXP:
  for (t = rval; t != R_NilValue; t = CDR(t), x = CDR(x)) {
    SETCAR(t, CAR(x));
    SET_TAG(t, TAG(x));
  }
  break;
case VECSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      SET_VECTOR_ELT(rval, i, VECTOR_ELT(x, i));
      if (xnames != R_NilValue)
        SET_STRING_ELT(names, i, STRING_ELT(xnames, i));
    }
    break;
case RAWSXP:
  for (i = 0; i < len; i++)
    if (i < lenx) {
      RAW(rval)[i] = RAW(x)[i];
      if (xnames != R_NilValue)
        SET_STRING_ELT(names, i, STRING_ELT(xnames, i));
    }
    else
      RAW(rval)[i] = (Rbyte) 0;
    break;
default:
  UNIMPLEMENTED_TYPE("length<-", x);
}
if (isVector(x) && xnames != R_NilValue)
  setAttrib(rval, R_NamesSymbol, names);
// *not* keeping "class": in line with  x[1:k]
UNPROTECT(2);
return rval;
  }
