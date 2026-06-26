#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXBUF 1000
#define MAXSEQ 5000

int main(int argc, char *argv[])
{
        FILE *fin;
        char buffer[MAXBUF];
        char bufseq[MAXSEQ];
        fpos_t pos;
        int aaa, aab;
	int output, upcut;

        if(argc!=2)
        {
                printf("takeseq.fas.x : a program to convert fasta format to takezaki format\n"); fflush(stdout);
                printf("use           : takeseq.fas.x fastafile\n"); fflush(stdout);
                exit(1);
        }

        if((fin=fopen(argv[1], "r"))==NULL)
        {
                printf("No fastafile\n"); fflush(stdout);
                exit(1);
        }

	while(fgets(buffer, MAXBUF, fin)!=NULL)
	{
		if(buffer[0]=='>')
		{
			fgetpos(fin, &pos);

			for(aaa=1; (buffer[aaa]!='\n' && buffer[aaa]!=' '); aaa++)
			{
				printf("%c", buffer[aaa]); fflush(stdout);
			}
			for(aaa; aaa<20; aaa++)
			{
				printf(" "); fflush(stdout);
			}

			aab=0;
			while(fgets(buffer, MAXBUF, fin)!=NULL)
			{
				if(buffer[0]!='>')
				{
					for(aaa=0; buffer[aaa]!='\n'; aaa++)
					{
						if(!isspace(buffer[aaa]))
						{
							bufseq[aab]=buffer[aaa];
							++aab;
						}
					}
					bufseq[aab]='\0';
				}
				else if(buffer[0]=='>')
				{
					goto owari1;
				}
			}

			owari1: ;
       		        fsetpos(fin, &pos);

			output=0;
			for(aaa=0; bufseq[aaa]!='\0'; aaa++)
			{
				if(isalpha(bufseq[aaa]))
				{
					printf("%c", bufseq[aaa]);
				}
			}
			printf("\n"); fflush(stdout);

		}
        }
        fclose(fin);
}
