#include "list.h"
#include <stdlib.h>
#include <stdio.h>

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
            //          printf("got list %li\n",generated_list.entries);
        } 
    }
    return generated_list;
}

void reallocateList(list_t *list,
        const size_t length) {
    entry_t *old_entries = list->entries;
    if (length>0) {
        if (list->length==0) 
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
    /*    if ((list->entries!=old_entries)&&(old_entries!=0)) {
          if (list->entries!=0)
          printf("moved %li by realloc to %li\n",old_entries,list->entries);
          else
          printf("freed %li by realloc\n", old_entries);
          }*/
}


void freeList(list_t *list) {
    if (list->entries!=0) {
        //        printf("deleting list %li\n",list->entries);
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
    else
        return 0;
}

void setEntry(const list_t list, const index_t index, const entry_t entry) {
    if ((index>=0)&&(index<getLength(list)))
        list.entries[index]=entry;
}


list_t initList(entry_t *entries, const size_t length) {
    index_t i;
    list_t l = allocateList(length);
    if (l.entries!=0) 
        for (i=0; i<(index_t)length; i++) 
            setEntry(l,i,entries[i]);
    return l;
}

list_t initListFromTo(const entry_t start, const entry_t stop) {
    index_t i;
    size_t length;
    entry_t c;
    int incr;
    if (stop>=start) {
        length=(size_t) (stop-start+1);
        incr=1;
    } else {
        length=(size_t) (start-stop+1);
        incr=-1;
    }
    list_t l = allocateList(length);
    if (l.entries!=0) 
        for (i=0,c=start; i<length; i++, c+=incr) 
            setEntry(l,i,c);
    return l;
}

list_t initConstantList(const entry_t c, const size_t length){
    entry_t i;
    list_t l = allocateList(length);
    if (l.entries!=0) 
        for (i=0; i<length; i++)
            setEntry(l,i,c);
    return l;
}


list_t duplicateList(const list_t list_in) {
    entry_t i;
    list_t list_out=emptyList();
    list_out = allocateList(getLength(list_in));
    for (i=0; i<getLength(list_out); i++) 
        setEntry(list_out,i,getEntry(list_in,i));
    return list_out;
}

list_t mergeLists(const list_t list1, const list_t list2) {
    list_t list_out;
    entry_t i,j;
    list_out = allocateList(getLength(list1)+ getLength(list2));
    if (list_out.entries!=0) {
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
        setEntry(new_list,i,getEntry(list,getEntry(indices,i)));
    }
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
    const index_t i=(index_t)getLength(*list);
    reallocateList(list,getLength(*list)+1);
    if (getLength(*list)>(size_t)i) {
        setEntry(*list,i,entry);
    }
}

void removeEntryFromList(list_t *list, const index_t index) {
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
        if (getEntry(*list,i)!=entry) 
            setEntry(*list, j++, getEntry(*list,i));
    }
    reallocateList(list,j);
}

void appendListToList(list_t *list1, const list_t list2) {
    index_t i,j;
    size_t siz_old = getLength(*list1);
    siz_old=list1->length;
    reallocateList(list1, getLength(*list1) + getLength(list2));
    if (getLength(*list1)>siz_old) {
        for (i=siz_old, j=0; i<list1->length; i++, j++) 
            setEntry(*list1,i,getEntry(list2,j));
    }
}

void removeEntryListFromList(list_t *list, const list_t indices) {
    index_t i,j;
    for (i=j=0; i<getLength(*list); i++)  
        if (notInList(i,indices)) 
            setEntry(*list, j++, getEntry(*list,i));
    reallocateList(list,j);
}

void removeValueListFromList(list_t *list, const list_t excl_list) {
    index_t i,j,k;
    int keep;
    for (j=i=0; i<getLength(*list); i++) {
        keep=1;
        for (k=0; k<getLength(excl_list); k++) 
            keep=(keep)&&(getEntry(*list,i)!=getEntry(excl_list,k));
        if (keep)
            setEntry(*list, j++, getEntry(*list,i));
    }
    reallocateList(list,j);
}

void reverseList(list_t * const list) {
    index_t i,j;
    entry_t v;
    const cnt = getLength(*list)/ 2;
    if (cnt>0)
        for (i=0, j=getLength(*list)-1; i<cnt; i++, j--) {
            v=getEntry(*list,i);
            setEntry(*list,i,getEntry(*list,j));
            setEntry(*list,j,v);
        }
}

int inList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) 
        if (getEntry(list,i)==entry)
            return 1;
    return 0;
}

entry_t findValueInList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) 
        if (entry==getEntry(list,i))
            return i;
    return getLength(list);  
}

void uniquefyListEntries(list_t *list)  {
    index_t i,j,k;
    int keep;
    k=0;
    for (j=0; j<getLength(*list); j++) {
        keep=1;
        for (i=0; (i<k)&&(keep); i++)  
            keep = (keep)&&(list->entries[j]!=list->entries[i]);
        if (keep) {
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
    for (i=0; i<getLength(value_list); i++) 
        appendToList(&l,findValueInList(
                    getEntry(value_list,i),list));
    return l;
}

int notInList(const entry_t entry, const list_t list) {
    index_t i;
    for (i=0; i<getLength(list); i++) 
        if (getEntry(list,i)==entry)
            return 0;
    return 1;
}

void printList(list_t const list) {
    entry_t i;
    printf("[list]_%d=[",list.length);
    if (getLength(list)>0)
        printf("%d",getEntry(list,0));
    for (i=1; i<list.length; i++) 
        printf(", %d",getEntry(list,i));
    printf("]\n");
}

