#include <stdio.h>
#include <stdlib.h>
#include "list.h"

int main(char **argv, int argc) {
   list_t l1=emptyList();
   list_t l2=emptyList();
   list_t l3=emptyList();
   int x[]={0, 2, 4, 6};

   printf("empty list:\n");
   printList(l1);
   freeList(&l1);

   printf("constant list with 10 1 entries:\n");
   l1=initConstantList(1,10);
   printList(l1);
   freeList(&l1);

   printf("list from 2 to 3:\n");
   l1=initListFromTo(2,3);
   printList(l1);
   printf("list from 7 to 1:\n");
   l2=initListFromTo(7,1);
   printList(l2);
   printf("duplicate list from 2 to 3:\n");
   l3=duplicateList(l1);
   printList(l3);
   freeList(&l3);

   printf("merge list from 1 to 7 with list from 2 to 3:\n");
   l3=mergeLists(l1,l2);
   printList(l3);

   printf("remove entry 6 list from list:\n");
   removeValueFromList(&l3,6);
   printList(l3);

   printf("remove entries [2, 3] from list:\n");
   removeValueListFromList(&l3,l1);
   printList(l3);

   printf("reverse list:\n");
   reverseList (&l3);
   printList(l3);

   printf("append entry 8 to list:\n");
   appendToList(&l3,8);
   printList(l3);

   printf("is 8 not in list?: %d\n",notInList(8,l3));
   printf("is 3 not in list?: %d\n",notInList(3,l3));

   printf("remove item 4 from list\n");
   removeEntryFromList(&l3,4);
   printList(l3);

   printf("remove item 1 from list\n");
   removeEntryFromList(&l3,1);
   printList(l3);

   printf("remove item 0 from list\n");
   removeEntryFromList(&l3,0);
   printList(l3);

   printf("...taking a longer list\n");
   printList(l2);
   freeList(&l3);
   l3=initList(x,4);
   printf("removing index list ");
   printList(l3);
   removeEntryListFromList(&l2,l3);
   printList(l2);
   freeList(&l1);
   freeList(&l2);
   freeList(&l3);

   l1=initListFromTo(1,3);
   l2=initListFromTo(0,5);
   printf("...taking a longer list\n");
   printList(l1);
   printf("finding indices of values in list");
   printList(l2);
   l3=findValueListInList(l2,l1);
   printList(l3);
   

   freeList(&l1);
   freeList(&l2);
   freeList(&l3);
   
}

// gcc list.c test_list.c && ./a.out
