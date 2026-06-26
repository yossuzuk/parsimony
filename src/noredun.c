#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXNAME 1000
#define MAXSEQ 10000
#define MAXENT 3000
#define MAXBUF 10000

main(int argc, char *argv[])
{
	struct entries {char name[MAXNAME]; char seq[MAXSEQ];};
	struct entries buf1;
	struct entries buf2;
	FILE *fin;
	char buffer1[MAXBUF], buffer2[MAXBUF];
	char same;
	int aa, ac, aaa, aab, aac, aad;
	int n; /* Variable for comparison. */
	int s=0; /* Variable for entry number. */
	int bn, bs; /* Variable for buffer structure. */
	int compare;
	int p, q;
	int output;
        fpos_t pos1, pos2;



	if(argc!=2) /* Error message for command. */
	{
		printf("nonredunfasta.x : a program for eliminating redundant sequences leaving the last one\n");
		printf("use             : nonredunfasta.x filename\n");
		exit(1);
	}



	buf1.name[0]='\0';
	buf1.seq[0]='\0';
	buf2.name[0]='\0';
	buf2.seq[0]='\0';



        if((fin=fopen(argv[1], "r"))==NULL)
        {
                printf("No fastafile\n"); fflush(stdout);
                exit(1);
        }
	while(fgets(buffer1, MAXBUF, fin)!=NULL)
	{
		if(buffer1[0]=='>')
		{
			fgetpos(fin, &pos1);
			aab=0;
			for(aaa=1; buffer1[aaa]!='\n'; aaa++)
			{
				buf1.name[aab]=buffer1[aaa];
				++aab;
			}
			buf1.name[aab]='\0';

			aab=0;
			while(fgets(buffer1, MAXBUF, fin)!=NULL)
			{
				if(buffer1[0]!='>')
				{
					for(aaa=0; !isspace(buffer1[aaa]); aaa++)
					{
						if(!isspace(buffer1[aaa]))
						{
							buf1.seq[aab]=buffer1[aaa];
							++aab;
						}
					}
					buf1.seq[aab]='\0';
				}
				else if(buffer1[0]=='>')
				{
					goto owari1;
				}
			}
			owari1: ;
       		        fsetpos(fin, &pos1);




			compare=0;
			while(fgets(buffer2, MAXBUF, fin)!=NULL)
			{
				if(buffer2[0]=='>')
				{
					fgetpos(fin, &pos2);
					aab=0;
					for(aaa=1; buffer2[aaa]!='\n'; aaa++)
					{
						buf2.name[aab]=buffer2[aaa];
						++aab;
					}
					buf2.name[aab]='\0';

					aab=0;
					while(fgets(buffer2, MAXBUF, fin)!=NULL)
					{
						if(buffer2[0]!='>')
						{
							for(aaa=0; !isspace(buffer2[aaa]); aaa++)
							{
								if(!isspace(buffer2[aaa]))
								{
									buf2.seq[aab]=buffer2[aaa];
									++aab;
								}
							}
							buf2.seq[aab]='\0';
						}
						else if(buffer2[0]=='>')
						{
							goto owari2;
						}
					}
					owari2: ;
					fsetpos(fin, &pos2);

					if(strcmp(buf1.seq, buf2.seq)==0)
					{
						compare=1;
						goto owari3;
					}
				}
			}

			owari3: ;
       			fsetpos(fin, &pos1);

			if(compare==0)
			{
				printf(">%s\n", buf1.name);

				output=0;
				for(aaa=0; buf1.seq[aaa]!='\0'; aaa++)
				{
					if(isalpha(buf1.seq[aaa]))
					{
						if(output>0 && output%60==0)
						{
							printf("\n"); fflush(stdout);
						}
						++output;
						printf("%c", buf1.seq[aaa]); fflush(stdout);
					}
				}
				printf("\n"); fflush(stdout);
			}
		}
	}
	fclose(fin);
}
