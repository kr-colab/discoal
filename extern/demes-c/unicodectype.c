/*
   This file and unicodetype_db.h were copied from the cpython 3.9.4 sources:
   https://github.com/python/cpython/blob/main/Objects/unicodectype.c
   https://github.com/python/cpython/blob/main/Objects/unicodetype_db.h
  
   These are liberally licensed:
   https://github.com/python/cpython/blob/main/LICENSE
  
   Only the necessary functionality to check for python identifiers
   has been retained (_PyUnicode_IsXidStart() and _PyUnicode_IsXidContinue).
*/

/*
   Unicode character type helpers.

   Written by Marc-Andre Lemburg (mal@lemburg.com).
   Modified for Python 2.0 by Fredrik Lundh (fredrik@pythonware.com)

   Copyright (c) Corporation for National Research Initiatives.

*/

#include <stdint.h>

typedef uint32_t Py_UCS4;

#define ALPHA_MASK 0x01
#define DECIMAL_MASK 0x02
#define DIGIT_MASK 0x04
#define LOWER_MASK 0x08
#define LINEBREAK_MASK 0x10
#define SPACE_MASK 0x20
#define TITLE_MASK 0x40
#define UPPER_MASK 0x80
#define XID_START_MASK 0x100
#define XID_CONTINUE_MASK 0x200
#define PRINTABLE_MASK 0x400
#define NUMERIC_MASK 0x800
#define CASE_IGNORABLE_MASK 0x1000
#define CASED_MASK 0x2000
#define EXTENDED_CASE_MASK 0x4000

typedef struct {
    /*
       These are either deltas to the character or offsets in
       _PyUnicode_ExtendedCase.
    */
    const int upper;
    const int lower;
    const int title;
    /* Note if more flag space is needed, decimal and digit could be unified. */
    const unsigned char decimal;
    const unsigned char digit;
    const unsigned short flags;
} _PyUnicode_TypeRecord;

#include "unicodetype_db.h"

static const _PyUnicode_TypeRecord *
gettyperecord(Py_UCS4 code)
{
    int index;

    if (code >= 0x110000)
        index = 0;
    else
    {
        index = index1[(code>>SHIFT)];
        index = index2[(index<<SHIFT)+(code&((1<<SHIFT)-1))];
    }

    return &_PyUnicode_TypeRecords[index];
}

/* Returns 1 for Unicode characters having the XID_Start property, 0
   otherwise. */

int _PyUnicode_IsXidStart(Py_UCS4 ch)
{
    const _PyUnicode_TypeRecord *ctype = gettyperecord(ch);

    return (ctype->flags & XID_START_MASK) != 0;
}

/* Returns 1 for Unicode characters having the XID_Continue property,
   0 otherwise. */

int _PyUnicode_IsXidContinue(Py_UCS4 ch)
{
    const _PyUnicode_TypeRecord *ctype = gettyperecord(ch);

    return (ctype->flags & XID_CONTINUE_MASK) != 0;
}
