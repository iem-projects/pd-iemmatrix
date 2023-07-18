#ifndef QHULL_ENTRY_H
#define QHULL_ENTRY_H

#include <sys/types.h>
#include <stdio.h>


/*
 *  variable entry types for
 *  list operations in zhull
 *
 * Copyright (c) 2012, IOhannes zmoelnig,
 * with friendly help from
 * IEM, Graz, Austria
 *
 *
 */



typedef size_t index_t;

typedef union {
  index_t i;
  void*p;
} entryvalu_t;

typedef enum {
  INDEX,
  POINTER,
  INVALID
} entrytype_t;

typedef struct entry_ {
  entrytype_t typ;
  entryvalu_t val;
} entry_t;

static
void entry_setIndex(entry_t*e, index_t i)
{
  e->typ=INDEX;
  e->val.i=i;
}
static
void entry_setPointer(entry_t*e, void*p)
{
  e->typ=POINTER;
  e->val.p=p;
}
static
entry_t entry_makeIndex(index_t i)
{
  entry_t result;
  entry_setIndex(&result, i);
  return result;
}
static
entry_t entry_makePointer(void*p)
{
  entry_t result;
  entry_setPointer(&result, p);
  return result;
}

static
index_t entry_getIndex(const entry_t*e)
{
  return (INDEX==e->typ)?e->val.i:0;
}
static
void*entry_getPointer(const entry_t*e)
{
  return (POINTER==e->typ)?e->val.p:0;
}
static
int entry_equals(const entry_t*e1, const entry_t*e2)
{
  if(e1->typ!=e2->typ) {
    return 0;
  }
  switch(e1->typ) {
  case INDEX:
    return (e1->val.i == e2->val.i);
  case POINTER:
    return (e1->val.p == e2->val.p);
  default:
    return 0;
  }
  return 0;
}

static
void print_entry(const entry_t e)
{
  switch(e.typ) {
  case INDEX:
    printf("%lu", (unsigned long)(e.val.i));
    return;
  case POINTER:
    printf("0x%p", e.val.p);
    return;
  default:
    printf("<unknown>");
    return;
  }
}





#endif
