#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

// Global variables for counting the number of calling quicksort
int quicksort_count = 0;

typedef struct item
{
    double dataPivot;
    double dataLeft;
    double dataRight;
}Item;

void swap(Item *a, Item *b) 
{
    Item temp = *a;
    *a = *b;
    *b = temp;
}

// Divide array into three parts and compare each element with pivot.
// After the comparison, the pivot would be put at the right position.
int partition(Item arr[], int low, int high) 
{
    double pivot = arr[high].dataPivot;
    int i = low - 1;

    for (int j = low; j <= high - 1; j++) 
    {
        if (arr[j].dataPivot < pivot) 
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return i + 1;
}

// Quicksort is conducted only when the high number is bigger than the low number.
void quicksort(Item arr[], int low, int high) 
{
    if (low < high) 
    {
        int pi = partition(arr, low, high);

        quicksort(arr, low, pi - 1);
        quicksort(arr, pi + 1, high);
        quicksort_count += 2;
    }
}


int main(int argc, char *argv[])
{
    FILE *file;
    char *filename;
    int num_elements = 0;
    int size;
    double left, pivot, right;

    filename = argv[1];

    // Open the file.
    file = fopen(filename, "r");

    // Read the number of items
    if (fscanf(file, "%d", &size) != 1) 
    {
        printf("Error reading size from the file.\n");
        fclose(file);
        return 1;
    }

    Item * itemArray = (Item *)malloc(size * sizeof(Item));


    // Read all items and put into struct array.
    while (fscanf(file, "%lf %lf %lf", &itemArray[num_elements].dataLeft, &itemArray[num_elements].dataPivot, &itemArray[num_elements].dataRight) != EOF) 
    {
        num_elements++;
    }
    fclose(file);

    // Conduct quicksort
    quicksort(itemArray, 0, size-1);
    quicksort_count += 1;



    // Print sorted array according "p" 
    if ((argc > 2) && (strcmp(argv[2], "p") == 0)) 
    {
        for (int i = 0; i < size; i++) 
        {
            printf("%.4f\t%.4f\t%.4f\n", itemArray[i].dataLeft, itemArray[i].dataPivot, itemArray[i].dataRight);
        }
        printf("%d\n", quicksort_count);
    }
    else
    {
        printf("%d\n", quicksort_count);
    }
    
    // Release the memory
    free(itemArray);

    return 0;

}