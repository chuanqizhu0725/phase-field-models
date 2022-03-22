#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i, j, l, k;
    int ND = 10;

    // int(*arr)[ND][ND][ND] = malloc(sizeof(*arr));

    // for (i = 0; i <= ND - 1; i++)
    // {
    //     for (j = 0; j <= ND - 1; j++)
    //     {
    //         for (l = 0; l <= ND - 1; l++)
    //         {
    //             (*arr)[i][j][l] = i * j * l;
    //         }
    //     }
    // }

    // int arr[ND][ND][ND];

    // for (i = 0; i <= ND - 1; i++)
    // {
    //     for (j = 0; j <= ND - 1; j++)
    //     {
    //         for (l = 0; l <= ND - 1; l++)
    //         {
    //             arr[i][j][l] = i * j * l;
    //         }
    //     }
    // }

    double(*arr)[3][ND][ND] = calloc(ND, sizeof *arr);
    for (k = 0; k <= 2; k++)
    {
        for (i = 0; i <= ND - 1; i++)
        {
            for (j = 0; j <= ND - 1; j++)
            {
                for (l = 0; l <= ND - 1; l++)
                {
                    arr[k][i][j][l] = 1000.0;
                }
            }
        }
    }
    free(arr);

    printf("%e\n", arr[2][ND - 1][ND - 1][ND - 1]);
    printf("hellw, c!\n");
    return 0;
}