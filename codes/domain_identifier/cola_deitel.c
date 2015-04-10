#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*
  struct queueNode
  {
  int data;
  struct queueNode *nextPtr;
  };
*/

int DisEmpty(struct queueNode *headPtr)
{
  return headPtr == NULL;
}


void Dencola(struct queueNode **headPtr,struct queueNode **tailPtr,int value)
{
  
  struct queueNode *newPtr=NULL;
  
  newPtr = (struct queueNode *) malloc(sizeof(struct queueNode));
  
  if(newPtr != NULL)
    {
      newPtr->data = value;
      newPtr->nextPtr = NULL;
      
      if(DisEmpty(*headPtr))
	*headPtr = newPtr;
      else
	(*tailPtr)->nextPtr = newPtr;
      
      *tailPtr = newPtr;
      
    }
  else
    printf("There is no memory available to storage data in tail\n");
  
}

int Ddesencola(struct queueNode **headPtr,struct queueNode **tailPtr)
{
  
  int value;
  struct queueNode *tempPtr=NULL;
  
  value = (*headPtr)->data;
  tempPtr = *headPtr;
  *headPtr= (*headPtr)->nextPtr;
  
  if(*headPtr == NULL)
    *tailPtr = NULL;
  
  free(tempPtr);
  
  return value;
  
}

void Dprintcola(struct queueNode *currentPtr)
{
  
  if(currentPtr == NULL)
    printf("queue is empty\n");
  else
    {
      printf("queue is in:\n");
      
      while(currentPtr != NULL)
	{
	  printf("%d --> ",currentPtr->data);
	  currentPtr = currentPtr->nextPtr;
	}
      
      printf("NULL\n\n");
    }
}


void Dsizecola(struct queueNode *currentPtr,int *count)
{
  
  while(currentPtr != NULL)
    {
      (*count) = (*count) + 1;
      currentPtr = currentPtr->nextPtr;
    }
  
}


/*
  int main()
  {
  
  int i,desenc,enc;
  struct queueNode *headPtr=NULL;
  struct queueNode *tailPtr=NULL;
  
  for(i=0; i<20; i++)
  {
  printf("encolando %d\n",i);
  encola(&headPtr,&tailPtr,i);
  }
  
  printf("********************\n");
  getchar();
  
  printcola(headPtr);
  
  
  for(i=0; i<10; i++)
  {
  desenc=desencola(&headPtr, &tailPtr);
  printf("Desencolando %d\n",desenc);
  }
  
  printf("********************\n");
  getchar();
  
  printcola(headPtr);
  
  for(i=100; i<110; i++)
  {
  printf("encolando %d\n",i);
  encola(&headPtr,&tailPtr,i);
  }
  
  printf("********************\n");
  getchar();
  
  printcola(headPtr);
  
  }
*/
