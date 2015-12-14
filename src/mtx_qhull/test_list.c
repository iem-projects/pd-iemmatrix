#include <stdio.h>
#include <stdlib.h>
#include "list.h"

int main(char **argv, int argc)
{
  list_t l1=emptyList();
  list_t l2=emptyList();
  list_t l3=emptyList();
  index_t x[]= {0, 2, 4, 6};
  printf("\nempty list:\n");
  printList(l1);
  freeList(&l1);

  printf("\nconstant list with 10 1 entries:\n");
  l1=initConstantList(entry_makeIndex(1),10);
  printList(l1);
  freeList(&l1);

  printf("\nlist from 2 to 3:\n");
  l1=initListFromTo(2,3);
  printList(l1);

  printf("\nlist from 7 to 1:\n");
  l2=initListFromTo(7,1);
  printList(l2);

  printf("\nduplicate list from 2 to 3:\n");
  l3=duplicateList(l1);
  printList(l3);
  freeList(&l3);

  printf("\nmerge list from 2..3 with list 7..1:\n");
  l3=mergeLists(l1,l2);
  printList(l3);

  printf("\nremove value 6 list from list:\n");
  removeValueFromList(&l3,entry_makeIndex(6));
  printList(l3);

  printf("\nremove values [2, 3] from list:\n");
  removeValueListFromList(&l3,l1);
  printList(l3);

  printf("\nreverse list:\n");
  reverseList (&l2);
  printList(l2);

  printf("\nappend entry 8 to list:\n");
  appendToList(&l3,entry_makeIndex(8));
  printList(l3);

  printf("\nis 8 not in list?: %d\n",notInList(entry_makeIndex(8),l3));
  printf("is 3 not in list?: %d\n",notInList(entry_makeIndex(3),l3));

  printf("\nremove index 4 from list\n");
  removeIndexFromList(&l3,4);
  printList(l3);

  printf("\nremove index 1 from list\n");
  removeIndexFromList(&l3,1);
  printList(l3);

  printf("\nremove index 0 from list\n");
  removeIndexFromList(&l3,0);
  printList(l3);

  printf("\n...taking a longer list\n");
  printList(l2);
  freeList(&l3);
  l3=initListIndex(x,4);
  printf("\nremoving index list ");
  printList(l3);
  removeEntryListFromList(&l2,l3);
  printList(l2);
  freeList(&l1);
  freeList(&l2);
  freeList(&l3);

  l1=initListFromTo(1,3);
  l2=initListFromTo(0,5);
  printf("\n...taking a longer list\n");
  printList(l1);
  printf("\nfinding indices of values in list");
  printList(l2);
  l3=findValueListInList(l2,l1);
  printList(l3);


  freeList(&l1);
  freeList(&l2);
  freeList(&l3);

}

// gcc list.c test_list.c && ./a.out
