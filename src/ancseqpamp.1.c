#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define MAXBUF 5000

int main(int argc, char *argv[])
{
	FILE *fin, *fin2, *fin3;
	FILE *fout, *fout2, *fout3;
	int aaa, aab, aac;
	int baa, bab, bac;
	int caa, cab, cac;
	int daa, dab, dac;
	int zaa, zab, zac;
	int oyanum;
	char buffer[MAXBUF];
	char buffer2[MAXBUF];
	char buffer3[MAXBUF];
	char oyaseq[MAXBUF];
	char oyanuma[20];
	char koname[20], seqname[20];

	if(argc!=4)
	{
		printf("ancseqpamp.x : a program to infer ancestral sequences using pamp\n"); fflush(stdout);
		printf("use          : ancseqpamp.x oyakofile [seqfile] refseq\n"); fflush(stdout);
		exit(1);
	}

/****** iterations ******/
/******
	for(zaa=1000001; zaa<=1018443; zaa++)
******/
	for(zaa=1000001; zaa<=1019948; zaa++)
	{
		fout=fopen("wseq", "w");
		fprintf(fout, "   3   1698\n");

		fout2=fopen("wtree", "w");
		fprintf(fout2, "   3   1\n");
		fprintf(fout2, "(");

		fin=fopen(argv[1], "r");
		while(fgets(buffer, MAXBUF, fin)!=NULL)
		{
			for(aaa=0; !isdigit(buffer[aaa]); aaa++) { ; }
			aab=0;
			for(aaa; isdigit(buffer[aaa]); aaa++)
			{
				oyanuma[aab]=buffer[aaa];
				++aab;
			}
			oyanuma[aab]='\0';
			oyanum=atoi(oyanuma);

			if(oyanum==zaa)
			{



				for(aaa; isspace(buffer[aaa]); aaa++) { ; }
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) { koname[aab]=buffer[aaa]; ++aab; }
				koname[aab]='\0';
/******
printf("%s\n", koname); fflush(stdout);
******/
				fin2=fopen("seqfile", "r");
				while(fgets(buffer2, MAXBUF, fin2)!=NULL)
				{
					bab=0;
					for(baa=0; !isspace(buffer2[baa]); baa++) { seqname[bab]=buffer2[baa]; ++bab; }
					seqname[bab]='\0';

					if(strcmp(seqname, koname)==0)
					{
						fprintf(fout, "%s", buffer2);
						for(caa=0; koname[caa]!='\0'; caa++)
						{
							fprintf(fout2, "%c", koname[caa]);
						}
						fprintf(fout2, ",");
					}
				}
				fclose(fin2);
			}
		}
		fin2=fopen("seqfile", "r");
		while(fgets(buffer2, MAXBUF, fin2)!=NULL)
		{
			bab=0;
			for(baa=0; !isspace(buffer2[baa]); baa++) { seqname[bab]=buffer2[baa]; ++bab; }
			seqname[bab]='\0';

			if(strcmp(seqname, argv[3])==0)
			{
				fprintf(fout, "%s", buffer2);
				for(caa=0; argv[3][caa]!='\0'; caa++)
				{
					fprintf(fout2, "%c", argv[3][caa]);
				}
			}
		}
		fprintf(fout2, ");\n");
		fclose(fin);
		fclose(fin2);
		fclose(fout);
		fclose(fout2);




		system("./pamp pamp.ctl");



		fin3=fopen("mp", "r");
		while(fgets(buffer3, MAXBUF, fin3)!=NULL)
		{
			if(buffer3[0]=='n' && buffer3[1]=='o' && buffer3[2]=='d' && buffer3[3]=='e')
			{
				for(daa=10; isspace(buffer3[daa]); daa++) { ; }

				baa=0;
				for(daa; (buffer3[daa]!='\n' && buffer3[daa]!='\0'); daa++)
				{
					if(isalpha(buffer3[daa]))
					{
						oyaseq[baa]=tolower(buffer3[daa]);
						++baa;
					}
					oyaseq[baa]='\0';
				}
			}
		}
		fclose(fin3);



		fout3=fopen("oya", "w");
		fprintf(fout3, "seq%d         ", zaa);
		for(daa=0; oyaseq[daa]!='\0'; daa++)
		{
			fprintf(fout3, "%c", oyaseq[daa]);
		}
		fprintf(fout3, "\n");

		fclose(fout3);



		system("cat seqfile oya > seqfile2");
		system("cp seqfile2 seqfile");






	}
}
