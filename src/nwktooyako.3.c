#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define MAXBUF 5000000
int g_buffer[5000000];

int main(int argc, char *argv[])
{
	FILE *fin;
	FILE *fout;
	static char buffer[MAXBUF];
	int aaa, aab, aac;
	int intnum; /****** number name of interior branch ******/
	int sentou;
	char name1[200], name2[200], blen1[200], blen2[200];

	if(argc!=2)
	{
		printf("nwkoyako.3 : a program to extract oyako from nwkfile sequence name very long to accession only\n"); fflush(stdout);
		printf("use        : nwkoyako.x nwkfile\n"); fflush(stdout);
		exit(1);
	}

	if((fin=fopen(argv[1], "r"))==NULL)
	{
		printf("No nwkfile\n"); fflush(stdout);
		exit(1);
	}
	fout=fopen("wtree", "w");
	fgets(buffer, MAXBUF, fin);


	intnum=1000000;

	for(aaa=0; (buffer[aaa]!='\n' && buffer[aaa]!='\0'); aaa++)
	{

		if(buffer[aaa]=='(' && isalnum(buffer[aaa+1]))
		{

			sentou=aaa;
			name1[0]=name2[0]=blen1[0]=blen2[0]='\0';

			aab=0;
			for(aaa=aaa+1; buffer[aaa]!=':'; aaa++) { name1[aab]=buffer[aaa]; ++aab; } name1[aab]='\0';
			aab=0;
			for(aaa=aaa+1; buffer[aaa]!=','; aaa++) { blen1[aab]=buffer[aaa]; ++aab; } blen1[aab]='\0';

			if(buffer[aaa+1]=='(')
			{
				goto owari2;
			}
			else if(isalnum(buffer[aaa+1]))
			{
				aab=0;
				for(aaa=aaa+1; buffer[aaa]!=':'; aaa++) { name2[aab]=buffer[aaa]; ++aab; } name2[aab]='\0';
				aab=0;
				for(aaa=aaa+1; (buffer[aaa]!=')' && buffer[aaa]!=','); aaa++) { blen2[aab]=buffer[aaa]; ++aab; } blen2[aab]='\0';

				if(buffer[aaa]==')' || buffer[aaa]==',')
				{
					++intnum;

					printf("seq%d", intnum); fflush(stdout);
					for(aac=0; aac<10; aac++) { printf(" "); fflush(stdout); }
					for(aac=0; name1[aac]!='\0'; aac++) { printf("%c", name1[aac]); fflush(stdout); }
					for(aac; aac<=200; aac++) { printf(" "); fflush(stdout); }
					for(aac=0; blen1[aac]!='\0'; aac++) { printf("%c", blen1[aac]); fflush(stdout); }
					printf("\n"); fflush(stdout);

					printf("seq%d", intnum); fflush(stdout);
					for(aac=0; aac<10; aac++) { printf(" "); fflush(stdout); }
					for(aac=0; name2[aac]!='\0'; aac++) { printf("%c", name2[aac]); fflush(stdout); }
					for(aac; aac<=200; aac++) { printf(" "); fflush(stdout); }
					for(aac=0; blen2[aac]!='\0'; aac++) { printf("%c", blen2[aac]); fflush(stdout); }
					printf("\n"); fflush(stdout);

/******
					system("egrep '(%s %s %s)' h3. > temp.seq", name1, name2);

					system("./paml paml.ctl");

					exit(1);
******/








					if(buffer[aaa]==')')
					{

						for(aac=0; aac<sentou; aac++)
						{
							fprintf(fout, "%c", buffer[aac]);
						}
						fprintf(fout, "seq%d", intnum);
						for(aac=aaa+1; (buffer[aac]!='\n' && buffer[aac]!='\0'); aac++)
						{
							fprintf(fout, "%c", buffer[aac]);
						}
						fclose(fout);
						fclose(fin);

						system("cp wtree rtree");

						fin=fopen("rtree", "r");
						fgets(buffer, MAXBUF, fin);
						fout=fopen("wtree", "w");

						aaa=-1;

						goto owari1;
					}
					else if(buffer[aaa]==',')
					{
						fclose(fout);
						fclose(fin);

						goto owari3;
					}

				}
			}
		}

		owari1: ;


		owari2: ;

	}

	owari3: ;

}
