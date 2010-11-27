#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
int main(void)
{
    FILE *fp;
    size_t count;
    const char *str = "hello\n";
 
    fp = fopen("sample.txt", "w");
    if(fp == NULL) {
        perror("failed to open sample.txt");
        return EXIT_FAILURE;
    }
    
    count = fwrite(str, 1, strlen(str), fp);
    count += fwrite(str, 1, strlen(str), fp);
	printf("Wrote %u bytes. fclose(fp) %s.\n", count, fclose(fp) == 0 ? "succeeded" : "failed");
    return EXIT_SUCCESS;
}

