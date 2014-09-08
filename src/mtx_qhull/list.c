#include "list.h"
#include <stdlib.h>
#include <stdio.h>

/*
 *  list operations for zhull
 *
 * Copyright (c) 2012, Franz Zotter,
 * with friendly help from
 * IOhannes zmoelnig
 * for variable entry types
 * in entry.h
 * IEM, Graz, Austria
 * 
 *
 */


// memory things:
list_t emptyList(void) {
    list_t generated_list;
    generated_list.length=0;
    generated_list.entries=0;
    return generated_list;
}

list_t allocateList(const size_t length) {
    list_t generated_list = emptyList();
    if (length>0) {
        generated_list.entries= (entry_t*) malloc(length*sizeof(entry_t));
        if (generated_list.entries!=0) {
            generated_list.length=length;
        } 
    }
    return generated_list;
}

void reallocateList(list_t *list,
        const size_t length) {
    entry_t *old_entries = list->entries;
    if (length>0) {
        if (getLength(*list)==0) 
            *list = allocateList(length);
        else {
            if (list->length != length)
                list->entries = (entry_t*) realloc(list->entries,length*sizeof(entry_t));
            if (list->entries!=0) 
                list->length=length;
            else 
                *list=emptyList();
        }
    } 
    else 
        freeList(list);
}


void freeList(list_t *list) {
    if (list->entries!=0) {
        free(list->entries);
    }
    list->entries=0;
    list->length=0;
}

// programming interface:

size_t getLength(const list_t list) {
    return list.length;
}

entry_t getEntry(const list_t list, const index_t index) {
    if ((index>=0)&&(index<getLength(list)))
        return list.entries[index];
    else {
      entry_t result={0};
        return result;
    }
}

void setEntry(const list_t list, const index_t index, 
        const entry_t entry) {
    if ((index>=0)&&(index<getLength(list)))
        list.entries[index]=entry;
}


list_t initList(const entry_t *entries, 
        const size_t length) {
    index_t i;
    list_t l = allocateList(length);
    if (getLength(l)!=0) 
        for (i=0; i<(index_t)length; i++) 
            setEntry(l,i,entries[i]);
    return l;
}

list_t initListIndex(const index_t *entries, const size_t length) {
  index_t i;
  list_t l = allocateList(length);
  if (getLength(l)!=0) 
    for (i=0; i<(index_t)length; i++) 
      setEntry(l,i,entry_makeIndex(entries[i]));
  return l;
}



list_t initListFromTo(const index_t start, const index_t stop) {
    entry_t e;
    index_t i;
    size_t length;
    index_t c;
    int incr;
    if (stop>=start) {
      length=(size_t) (stop-start+1);
        incr=1;
    } else {
      length=(size_t) (start-stop+1);
      incr=-1;
    }
    list_t l = allocateList(length);
    if (getLength(l)!=0) {
      for (i=0,c=start; i<length; i++, c+=incr) {
        entry_setIndex(&e, c);
        setEntry(l,i,e);
      }
    }
    return l;
}

list_t initConstantList(const entry_t c, const size_t length){
    index_t i;
    list_t l = allocateList(length);
    if (getLength(l)!=0) 
        for (i=0; i<length; i++)
            setEntry(l,i,c);
    return l;
}


list_t duplicateList(const list_t list_in) {
    index_t i;
    list_t list_out=emptyList();
    list_out = allocateList(getLength(list_in));
    for (i=0; i<getLength(list_out); i++) 
        setEntry(list_out,i,getEntry(list_in,i));
    return list_out;
}

list_t mergeLists(const list_t list1, const list_t list2) {
    list_t list_out;
    index_t i,j;
    list_out = allocateList(getLength(list1)+ getLength(list2));
    if (getLength(list_out)>=getLength(list1)) {
        for (i=0; i<getLength(list1); i++) 
            setEntry(list_out,i,getEntry(list1,i));
        for (j=0; i<getLength(list_out); i++, j++) 
            setEntry(list_out,i,getEntry(list2,j));
    }
    return list_out;
}

list_t getSubList(const list_t list, const list_t indices) {
    index_t i;
    list_t new_list = allocateList(getLength(indices));
    for (i=0; i<getLength(new_list); i++) {
        entry_t e1=getEntry(indices,i);
        setEntry(new_list,i,
                 getEntry(list,entry_getIndex(&e1))
                 );
    }
    return new_list;
}

list_t getSubListFromTo(const list_t list, const index_t start, 
        const index_t stop)  {
    list_t new_list=emptyList();
    index_t i,j;
    int incr;
    if ((start>0)&&(stop>0)&&(start<getLength(list))&&(stop<getLength(list))) {
        if (start>stop) {
            incr=-1;
            new_list=allocateList(start-stop+1);
        } else {
            incr=1;
            new_list=allocateList(start-stop+1);
        }
        for (j=start,i=0; i<getLength(new_list); i++, j+=incr) {
            setEntry(new_list,i,getEntry(list,j));
        }
    }
    return new_list;
}

void appendToList(list_t *list, const entry_t entry) {
    const size_t i=getLength(*list);
    reallocateList(list,getLength(*list)+1);
    if (getLength(*list)>i) {
        setEntry(*list,i,entry);
    }
}

void removeIndexFromList(list_t *list, const index_t index) {
    index_t i,j;
    for (i=j=0; i<getLength(*list); i++) {
        if (i!=index) 
            setEntry(*list,j++,getEntry(*list,i));
    }
    reallocateList(list,j);
}


void removeValueFromList(list_t *list, const entry_t entry) {
    index_t i,j;
    for (j=i=0; i<getLength(*list); i++) {
      entry_t e1=getEntry(*list,i);
      if (!entry_equals(&e1, &entry)) 
        setEntry(*list, j++, getEntry(*list,i));
    }
    reallocateList(list,j);
}

void appendListToList(list_t *list1, const list_t list2) {
    index_t i,j;
    const size_t siz_old = getLength(*list1);
    reallocateList(list1, getLength(*list1) + getLength(list2));
    for (i=siz_old, j=0; i<getLength(*list1); i++, j++)  
        setEntry(*list1,i,getEntry(list2,j));
}

void removeEntryListFromList(list_t *list, const list_t indices) {
    index_t i,j;
    for (i=j=0; i<getLength(*list); i++) {
      entry_t e={i};
      if (notInList(e,indices)) 
        setEntry(*list, j++, getEntry(*list,i));
    }
    reallocateList(list,j);
}

void removeValueListFromList(list_t *list, const list_t excl_list) {
    index_t i,j,k;
    int keep;
    for (j=i=0; i<getLength(*list); i++) {
        keep=1;
        for (k=0; k<getLength(excl_list); k++) {
          entry_t e1=getEntry(*list, i);
          entry_t e2=getEntry( excl_list, k);
          keep=(keep)&&(!entry_equals(&e1, &e2));
        }
        if (keep)
            setEntry(*list, j++, getEntry(*list,i));
    }
    reallocateList(list,j);
}

void reverseList(list_t * const list) {
    index_t i,j;
    entry_t v;
    const size_t cnt = getLength(*list)/ 2;
    if (cnt>0)
        for (i=0, j=getLength(*list)-1; i<cnt; i++, j--) {
            v=getEntry(*list,i);
            setEntry(*list,i,getEntry(*list,j));
            setEntry(*list,j,v);
        }
}

int inList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) {
      entry_t e1=getEntry(list,i);
      if(entry_equals(&e1, &entry))
        return 1;
    }
    return 0;
}

index_t findValueInList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) {
      entry_t e1=getEntry(list,i);
      if(entry_equals(&e1, &entry))
        return i;
    }
    return i;
}

void uniquefyListEntries(list_t *list)  {
    index_t i,j,k;
    k=0;
    for (j=0; j<getLength(*list); j++) {
        for (i=0; i<k; i++) {
          entry_t e1=list->entries[j];
          entry_t e2=list->entries[i];
          if(entry_equals(&e1, &e2))
            break;
        }
        if (i==k) {
            list->entries[i++]=list->entries[j];
            k++;
        }
    }
    reallocateList(list, k);
}

list_t findValueListInList(const list_t value_list,
        const list_t list) {
    list_t l=emptyList();
    index_t i,j;
    for (i=0; i<getLength(value_list); i++) {
      index_t idx=findValueInList(getEntry(value_list,i),list);
      entry_t e;
      entry_setIndex(&e, idx);
      appendToList(&l, e);
    }
    return l;
}

int notInList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) {
      entry_t e=getEntry(list, i);
      if (entry_equals(&e, &entry))
        return 0;
    }
    return 1;
}

void printList(list_t const list) {
    index_t i;
    const size_t len=getLength(list);
    printf("[list]_%d=[",(int)len);
    if (len>0) {
      print_entry(getEntry(list,0));
    } 
    for (i=1; i<len; i++) {
      printf(", ");
      print_entry(getEntry(list,i));
    }
    printf("]\n");
}

