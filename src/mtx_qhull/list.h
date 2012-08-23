#ifndef QHULL_LIST_H
#define QHULL_LIST_H

#include <sys/types.h>

typedef long int entry_t;
typedef long int index_t;

typedef struct list_ {
    entry_t *entries;
    size_t length;
} list_t;

// memory things:
list_t emptyList(void);
void freeList(list_t *list);

// programming interface:
size_t getLength(const list_t list);
entry_t getEntry(const list_t list, const index_t index);
void setEntry(const list_t list, const index_t index, const entry_t entry);
list_t initList(entry_t *entries, const size_t length);
list_t initListFromTo(const entry_t start, const entry_t stop);
list_t initConstantList(const entry_t c, const size_t length);
list_t duplicateList(const list_t list_in);
list_t mergeLists(const list_t list1, const list_t list2);
list_t getSubList(const list_t list, const list_t indices);
list_t getSubListFromTo(const list_t list, const index_t start, 
      const index_t stop);
void appendToList(list_t *list, const entry_t entry);
void removeValueFromList(list_t *list, const entry_t entry);
void removeEntryFromList(list_t *list, const index_t index);
void appendListToList(list_t *list1, const list_t list2);
void removeValueListFromList(list_t *list, const list_t excl_list);
void removeEntryListFromList(list_t *list, const list_t indices);
void reverseList(list_t * const list);
int inList(const entry_t entry, const list_t list);
int notInList(const entry_t entry, const list_t list);
list_t findValueListInList(const list_t value_list, const list_t list);
entry_t findValueInList(const entry_t entry, const list_t list);
void printList(const list_t list);

#endif /* QHULL_LIST_H */
